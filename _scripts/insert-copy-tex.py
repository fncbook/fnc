#!/usr/bin/env python3
import sys
from pathlib import Path

TAG = '<script defer src="https://cdn.jsdelivr.net/npm/katex@0/dist/contrib/copy-tex.min.js" crossorigin="anonymous"></script>\n'

build_dir = Path("_build/html")
files = list(build_dir.rglob("*.html"))

for path in files:
    text = path.read_text(encoding="utf-8")
    if "</head>" not in text:
        continue
    if "copy-tex" in text:
        continue  # idempotent
    path.write_text(text.replace("</head>", TAG + "</head>", 1), encoding="utf-8")

print(f"Injected copy-tex into {len(files)} files.")