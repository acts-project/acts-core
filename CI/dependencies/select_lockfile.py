#!/usr/bin/env python3

import os
import argparse
import json
import urllib.request
import urllib.error
import re
import subprocess
import hashlib
import tempfile
from pathlib import Path
from typing import Tuple, Dict, Optional
import contextlib

# Modify the default cache dir to use a temporary directory
DEFAULT_CACHE_SIZE_LIMIT = 1 * 1024 * 1024  # 1MB


def compute_cache_key(url: str) -> str:
    """Compute a cache key for a URL"""
    return hashlib.sha256(url.encode()).hexdigest()


def compute_cache_digest(cache_dir: Path) -> str:
    """Compute a digest of all cache files except digest.txt"""
    files = sorted(
        f
        for f in os.listdir(cache_dir)
        if (cache_dir / f).is_file() and f != "digest.txt"
    )

    digest = hashlib.sha256()
    for fname in files:
        fpath = cache_dir / fname
        digest.update(fname.encode())
        digest.update(str(fpath.stat().st_size).encode())
        digest.update(fpath.read_bytes())
    return digest.hexdigest()


def update_cache_digest(cache_dir: Path):
    """Update the cache digest file"""
    digest = compute_cache_digest(cache_dir)
    (cache_dir / "digest.txt").write_text(digest)


def prune_cache(cache_dir: Optional[Path], size_limit: int):
    """Prune the cache to keep it under the size limit"""
    if cache_dir is None or not cache_dir.exists():
        return

    # Get all cache files with their modification times
    cache_files = [
        (cache_dir / f, (cache_dir / f).stat().st_mtime)
        for f in os.listdir(cache_dir)
        if (cache_dir / f).is_file()
        and f != "digest.txt"  # Exclude digest from pruning
    ]
    total_size = sum(f.stat().st_size for f, _ in cache_files)

    if total_size <= size_limit:
        return

    # Sort by modification time (oldest first)
    cache_files.sort(key=lambda x: x[1])

    # Remove files until we're under the limit
    for file_path, _ in cache_files:
        if total_size <= size_limit:
            break
        total_size -= file_path.stat().st_size
        file_path.unlink()

    # Update digest after pruning
    update_cache_digest(cache_dir)


def fetch_github(base_url: str, cache_dir: Optional[Path], cache_limit: int) -> bytes:
    headers = {}
    token = os.environ.get("GITHUB_TOKEN")
    if token is not None and token != "":
        headers["Authorization"] = f"token {token}"

    with contextlib.ExitStack() as stack:
        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            cache_dir = Path(stack.enter_context(tempfile.TemporaryDirectory()))

        # Check cache first
        cache_key = compute_cache_key(base_url)
        cache_file = cache_dir / cache_key

        if cache_file.exists():
            print("Cache hit on", base_url)
            return cache_file.read_bytes()
        else:
            print("Cache miss on", base_url)

        try:
            req = urllib.request.Request(base_url, headers=headers)
            with urllib.request.urlopen(req) as response:
                content = response.read()

                # Write to cache
                cache_file.write_bytes(content)

                # Update digest after adding new file
                update_cache_digest(cache_dir)

                # Prune cache if necessary (this will update digest again if pruning occurs)
                prune_cache(cache_dir, cache_limit)

                return content
        except urllib.error.URLError as e:
            print(f"Failed to fetch from {base_url}: {e}")
            exit(1)
        except json.JSONDecodeError as e:
            print(f"Failed to parse JSON response: {e}")
            exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", type=str, required=True, help="Tag to use")
    parser.add_argument("--arch", type=str, required=True, help="Architecture to use")
    parser.add_argument(
        "--compiler-binary",
        type=str,
        default=os.environ.get("CXX"),
        help="Compiler to use (defaults to CXX environment variable if set)",
    )
    parser.add_argument(
        "--compiler",
        type=str,
        default=None,
        help="Compiler to use (defaults to compiler binary if set)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file to write lockfile to",
    )
    parser.add_argument(
        "--cache-dir",
        type=lambda x: Path(x).expanduser() if x else None,
        default=os.environ.get("LOCKFILE_CACHE_DIR"),
        help="Directory to use for caching (defaults to LOCKFILE_CACHE_DIR env var)",
    )
    parser.add_argument(
        "--cache-limit",
        type=int,
        default=int(os.environ.get("LOCKFILE_CACHE_LIMIT", DEFAULT_CACHE_SIZE_LIMIT)),
        help="Cache size limit in bytes (defaults to LOCKFILE_CACHE_LIMIT env var)",
    )
    args = parser.parse_args()

    print("Fetching lockfiles for tag:", args.tag)
    print("Architecture:", args.arch)

    base_url = f"https://api.github.com/repos/acts-project/ci-dependencies/releases/tags/{args.tag}"

    data = json.loads(fetch_github(base_url, args.cache_dir, args.cache_limit))

    lockfiles = parse_assets(data)

    print("Available lockfiles:")
    for arch, compilers in lockfiles.items():
        print(f"> {arch}:")
        for c, (n, _) in compilers.items():
            print(f"  - {c}: {n}")

    if args.arch not in lockfiles:
        print(f"No lockfile found for architecture {args.arch}")
        exit(1)

    if args.compiler_binary is not None:
        compiler = determine_compiler_version(args.compiler_binary)
        print("Compiler:", args.compiler_binary, f"{compiler}")
    elif args.compiler is not None:
        if not re.match(r"^([\w-]+)@(\d+\.\d+\.\d+)$", args.compiler):
            print(f"Invalid compiler format: {args.compiler}")
            exit(1)
        compiler = args.compiler
        print("Compiler:", f"{compiler}")
    else:
        compiler = None

    lockfile = select_lockfile(lockfiles, args.arch, compiler)

    print("Selected lockfile:", lockfile)

    if args.output:
        with open(args.output, "wb") as f:
            f.write(fetch_github(lockfile, args.cache_dir, args.cache_limit))


def parse_assets(data: Dict) -> Dict[str, Dict[str, Tuple[str, str]]]:
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]] = {}

    for asset in data["assets"]:
        url = asset["browser_download_url"]

        name = asset["name"]
        if not name.endswith(".lock") or not name.startswith("spack_"):
            continue

        m = re.match(r"spack_(.*(?:aarch64|x86_64))(?:_(.*))?\.lock", name)
        if m is None:
            continue

        arch, compiler = m.groups()
        compiler = compiler if compiler else "default"
        lockfiles.setdefault(arch, {})[compiler] = (name, url)

    return lockfiles


def select_lockfile(
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]], arch: str, compiler: Optional[str]
):
    # Default to the default lockfile
    _, lockfile = lockfiles[arch]["default"]

    if compiler is None:
        return lockfile

    # Extract compiler family and version
    compiler_family = compiler.split("@")[0]

    # Find all matching compiler families
    matching_compilers = {
        comp: ver
        for comp, ver in lockfiles[arch].items()
        if comp != "default" and comp.split("@")[0] == compiler_family
    }

    if matching_compilers:
        if compiler in matching_compilers:
            # Exact match found
            _, lockfile = matching_compilers[compiler]
        else:
            # Find highest version of same compiler family
            highest_version = max(
                matching_compilers.keys(),
                key=lambda x: [int(v) for v in x.split("@")[1].split(".")],
            )
            _, lockfile = matching_compilers[highest_version]

    return lockfile


def determine_compiler_version(binary: str):
    try:
        result = subprocess.run([binary, "--version"], capture_output=True, text=True)

        line = result.stdout.split("\n", 1)[0]
        print(line)
        if "clang" in line:
            compiler = "clang"
            if "Apple" in line:
                compiler = "apple-clang"
        elif "gcc" in line or "GCC" in line or "g++" in line:
            compiler = "gcc"
        else:
            print(f"Unknown compiler: {binary}")
            exit(1)

        m = re.search(r"(\d+\.\d+\.\d+)", line)
        if m is None:
            print(f"Failed to determine version for compiler: {binary}")
            exit(1)
        (version,) = m.groups()
        return f"{compiler}@{version}"

    except (subprocess.SubprocessError, FileNotFoundError):
        print(f"Failed to determine version for compiler: {binary}")
        exit(1)


if __name__ == "__main__":
    main()
