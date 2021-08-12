// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/detail/Aligned.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/ArenaCleaner.hpp"

// System include(s)
#include <cstddef>
#include <utility>
#include <map>
#include <shared_mutex>
#include <thread>

namespace Acts {
namespace Nmm {

struct CudaDeviceId {
		using valueType = int;

		explicit constexpr CudaDeviceId(valueType id) noexcept : id_{id} {}

		constexpr valueType value() const noexcept { return id_; }

	private:
		valueType id_;
};

namespace MemoryResource {

class DeviceMemoryResource {
	public:
		virtual ~DeviceMemoryResource() = default;

		void* allocate(std::size_t bytes, CudaStreamView stream = CudaStreamView{}) {
			return doAllocate(Nmm::detail::align_up(bytes, 8), stream);
		}

		void deallocate(void* p, std::size_t bytes, CudaStreamView stream = CudaStreamView{}) {
			doDeallocate(p, Nmm::detail::align_up(bytes, 8), stream);
		}

		bool isEqual(DeviceMemoryResource const& other) const noexcept {
			return doIsEqual(other);
		}

		virtual bool supportsStreams() const noexcept = 0;

		virtual bool supportsGetMemInfo() const noexcept = 0;

		std::pair<std::size_t, std::size_t> getMemInfo(CudaStreamView stream) const {
			return doGetMemInfo(stream);
		}

	private:
		virtual void* doAllocate(std::size_t bytes, CudaStreamView stream) = 0;

		virtual void doDeallocate(void* p, std::size_t bytes, CudaStreamView stream) = 0;

		virtual bool doIsEqual(DeviceMemoryResource const& other) const noexcept {
			return this == &other;
		}

		virtual std::pair<std::size_t, std::size_t> doGetMemInfo(CudaStreamView stream) const = 0;

};// class DeviceMemoryResource

class CudaMemoryResource final : public DeviceMemoryResource {
	public:
		CudaMemoryResource() = default;
		~CudaMemoryResource() = default;
		CudaMemoryResource(CudaMemoryResource const&) = default;
		CudaMemoryResource(CudaMemoryResource&&) = default;
		CudaMemoryResource& operator=(CudaMemoryResource const&) = default;
		CudaMemoryResource& operator=(CudaMemoryResource&&) = default;

		bool supportsStreams() const noexcept override { return false; }

		bool supportsGetMemInfo() const noexcept override { return true; }

	private:
		void* doAllocate(std::size_t bytes, CudaStreamView) override {
			void* p{nullptr};
			cudaMalloc(&p, bytes);
			return p;
		}

		void doDeallocate(void* p, std::size_t, CudaStreamView) override {
			cudaFree(p);
		}

		bool doIsEqual(DeviceMemoryResource const& other) const noexcept override {
			return dynamic_cast<CudaMemoryResource const*>(&other) != nullptr;
		}

		std::pair<size_t, size_t> doGetMemInfo(CudaStreamView) const override {
			std::size_t freeSize;
			std::size_t totalSize;

			cudaMemGetInfo(&freeSize, &totalSize);
			return std::make_pair(freeSize, totalSize);
		}

};// class CudaMemoryResource

template <typename Upstream>
class ArenaMemoryResource final : public DeviceMemoryResource {
	public:
		explicit ArenaMemoryResource(Upstream* upstreamMemoryResource, std::size_t initialSize = globalArena::defaultInitialSize, std::size_t maximumSize = globalArena::defaultMaximumSize)
			: globalArena_{upstreamMemoryResource, initialSize, maximumSize}
			{}

		ArenaMemoryResource(ArenaMemoryResource const&) = delete;
		ArenaMemoryResource& operator=(ArenaMemoryResource const&) = delete;

		bool supportsStreams() const noexcept override { return true; }

		bool supportsGetMemInfo() const noexcept override { return false; }

	private:
		using globalArena = detail::Arena::GlobalArena<Upstream>;
		using arena = detail::Arena::Arena<Upstream>;
		using readLock = std::shared_lock<std::shared_timed_mutex>;
		using writeLock = std::lock_guard<std::shared_timed_mutex>;

		void* doAllocate(std::size_t bytes, CudaStreamView stream) override {
			if(bytes <= 0) return nullptr;

			bytes = detail::Arena::alignUp(bytes);
			return getArena(stream).allocate(bytes);
		}

		void doDeallocate(void* p, std::size_t bytes, CudaStreamView stream) override {
			if(p == nullptr || bytes <= 0) return;

			bytes = detail::Arena::alignUp(bytes);
			if(!getArena(stream).deallocate(p, bytes, stream)) {
				deallocateFromOtherArena(p, bytes, stream);
			}
		}

		void deallocateFromOtherArena(void* p, std::size_t bytes, CudaStreamView stream) {
			stream.synchronize();

			readLock lock(mtx_);

			if(usePerThreadArena(stream)) {
				auto const id = std::this_thread::get_id();
				for(auto& kv : threadArenas_) {
					if(kv.first != id && kv.second->deallocate(p, bytes)) return;
				}
			} else {
				for(auto& kv : streamArenas_){
					if(stream != kv.first && kv.second.deallocate(p, bytes)) return;
				}
			}
			globalArena_.deallocate({p, bytes});
		}

		arena& getArena(CudaStreamView stream) {
			if(usePerThreadArena(stream)){
				return getThreadArena();
			} else {
				return getStreamArena(stream);
			}
		}

		arena& getThreadArena() {
			auto const id = std::this_thread::get_id();
			
			readLock lockRead(mtx_);
			auto const it = threadArenas_.find(id);
			if(it != threadArenas_.end()) {
				return *it->second;
			}

			writeLock lockWrite(mtx_);
			auto a = std::make_shared<arena>(globalArena_);
			threadArenas_.emplace(id, a);
			thread_local detail::Arena::ArenaCleaner<Upstream> cleaner{a};
			return *a;
		}

		arena& getStreamArena(CudaStreamView stream) {
			//if(usePerThreadArena(stream)){

				readLock lockRead(mtx_);
				auto const it = streamArenas_.find(stream.value());
				if (it != streamArenas_.end()) { return it->second; }

				writeLock lockWrite(mtx_);
				streamArenas_.emplace(stream.value(), globalArena_);
				return streamArenas_.at(stream.value());
			//}
		}

		std::pair<std::size_t, std::size_t> doGetMemInfo(CudaStreamView stream) const override {
			return std::make_pair(0, 0);
		}

		static bool usePerThreadArena(CudaStreamView stream) {
			return stream.is_per_thread_default();
		}

		globalArena globalArena_;
		std::map<std::thread::id, std::shared_ptr<arena>> threadArenas_;
		std::map<cudaStream_t, arena> streamArenas_;
		mutable std::shared_timed_mutex mtx_;
};// class ArenaMemoryResource

namespace detail {

inline DeviceMemoryResource* initialResource() {
	static CudaMemoryResource memoryResource{};
	return &memoryResource;
}

inline std::mutex& mapLock() {
	static std::mutex mapLockReturn;
	return mapLockReturn;
}

inline auto& getMap() {
	static std::map<CudaDeviceId::valueType, DeviceMemoryResource*> deviceIdToResource;
	return deviceIdToResource;
}

inline CudaDeviceId currentDevice() {
	int deviceId;
	cudaGetDevice(&deviceId);
	return CudaDeviceId{deviceId};
}

inline DeviceMemoryResource* getPerDeviceResource(CudaDeviceId id) {
	std::lock_guard<std::mutex> lock{detail::mapLock()};
	auto& map = detail::getMap();

	auto const found = map.find(id.value());

	return (found == map.end()) ? (map[id.value()] = detail::initialResource()) : found->second;
}

inline DeviceMemoryResource* setPerDeviceResource(CudaDeviceId id, DeviceMemoryResource* newMemoryResource) {
	std::lock_guard<std::mutex> lock{detail::mapLock()};
	auto& map = detail::getMap();
	auto const oldItr = map.find(id.value());

	auto oldMemoryResource = (oldItr == map.end()) ? detail::initialResource() : oldItr->second;
	map[id.value()] = (newMemoryResource == nullptr) ? detail::initialResource() : newMemoryResource;
	return oldMemoryResource;
}

inline DeviceMemoryResource* getCurrentDeviceResource() {
	return getPerDeviceResource(detail::currentDevice());
}

inline DeviceMemoryResource* setCurrentDeviceResource(DeviceMemoryResource* newMemoryResource) {
	return setPerDeviceResource(detail::currentDevice(), newMemoryResource);
}

}
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts