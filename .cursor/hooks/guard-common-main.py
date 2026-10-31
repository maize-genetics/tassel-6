#!/usr/bin/env python3
"""Deny edits that add JVM-only imports or references to commonMain/commonTest Kotlin files."""

import json
import re
import sys

FORBIDDEN = re.compile(
    r"^\s*import\s+(java|javax|kotlin\.jvm|htsjdk|com\.google\.common|org\.apache)\.",
    re.MULTILINE,
)
QUALIFIED = re.compile(r"\b(java|javax)\.(io|nio|util|lang|text|math|net)\.[A-Za-z_][\w.]*")


def allow():
    print(json.dumps({"permission": "allow"}))
    sys.exit(0)


try:
    event = json.load(sys.stdin)
except Exception:
    allow()

tool_input = event.get("tool_input") or {}
path = tool_input.get("path") or tool_input.get("file_path") or ""
if not re.search(r"/src/common(Main|Test)/", path) or not path.endswith(".kt"):
    allow()

texts = [tool_input.get(key) for key in ("contents", "new_string", "content")]
texts += [edit.get("new_string") for edit in tool_input.get("edits") or [] if isinstance(edit, dict)]
text = "\n".join(t for t in texts if isinstance(t, str))

hits = sorted({line.strip() for line in text.splitlines() if FORBIDDEN.match(line)} |
              {m.group(0) for m in QUALIFIED.finditer(text)})
if not hits:
    allow()

print(json.dumps({
    "permission": "deny",
    "user_message": "Blocked a JVM-only reference in commonMain: " + ", ".join(hits),
    "agent_message": (
        "commonMain must stay pure Kotlin (AGENTS.md hard constraint 1, ADR-0001). "
        "Found: " + ", ".join(hits) + ". Use kotlin.*, kotlinx.io, or kotlinx.coroutines, "
        "or move platform code to jvmMain behind expect/actual or an interface."
    ),
}))
sys.exit(0)
