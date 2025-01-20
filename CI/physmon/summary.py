#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import os
import csv

HERALD_URL = "https://acts-herald.app.cern.ch/view/{repo}/runs/{run_id}/artifacts/{artifact_name}/{path}"
IS_CI = "GITHUB_ACTIONS" in os.environ


parser = argparse.ArgumentParser()
parser.add_argument("results")
parser.add_argument("--html")
parser.add_argument("--md")
args = parser.parse_args()

re_title = re.compile(r'<p class="title">\s*(.*)\s*<\/p>', re.RegexFlag.MULTILINE)
re_check = re.compile(r'<a.*title="(.*)">\s*(.)\s*<\/a>', re.RegexFlag.MULTILINE)

summary = []

with open(args.results) as f:
    reader = csv.reader(f)
    for title, html_path, ec in reader:
        summary.append(
            {
                "title": title,
                "total": ec == "0",
                "path": html_path,
            }
        )

if args.html:
    with open(args.html, mode="w", encoding="utf-8") as f:
        f.write(
            """<!DOCTYPE html>
<html>
<head>
  <title>physmon summary</title>
  <meta charset="UTF-8">
</head>
<body>
  <h1>physmon summary</h1>
  <ul>
            """
        )

        for s in summary:
            if s["title"].startswith("Comparison"):
                f.write(
                    f"""
        <li>🔵 <a href="{s["path"]}">{s["title"]}</a></li>"""
                )
            else:
                f.write(
                    f"""
        <li>{"✅" if s["total"] else "🔴"} <a href="{s["path"]}">{s["title"]}</a></li>"""
                )

        f.write(
            """
      </ul>
    </body>
    </html>
            """
        )

if args.md:
    with open(args.md, mode="w", encoding="utf-8") as f:
        f.write("# physmon summary\n")
        for s in summary:
            if IS_CI:
                url = HERALD_URL.format(
                    repo=os.environ["GITHUB_REPOSITORY"],
                    run_id=os.environ["GITHUB_RUN_ID"],
                    artifact_name="physmon",
                    path=s["path"],
                )
            else:
                url = s["path"]

            if s["title"].startswith("Comparison"):
                f.write(f"  - 🔵️ [{s['title']}]({url})\n")
            else:
                f.write(f"  - {'✅' if s['total'] else '🔴'} [{s['title']}]({url})\n")
