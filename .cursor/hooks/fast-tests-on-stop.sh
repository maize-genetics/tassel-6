#!/usr/bin/env bash
# On agent stop, run the fast JVM tests if Kotlin sources changed, and ask the
# agent to fix failures. Does nothing until the tassel-core module exists.

input=$(cat)
status=$(printf '%s' "$input" | jq -r '.status // empty' 2>/dev/null)

if [[ "$status" != "completed" || ! -d tassel-core || ! -x ./gradlew ]]; then
  echo '{}'
  exit 0
fi

changed=$(git status --porcelain -- '*.kt' '*.kts' 2>/dev/null)
if [[ -z "$changed" ]]; then
  echo '{}'
  exit 0
fi

log=$(mktemp -t tassel-fast-tests)
if ./gradlew :tassel-core:jvmTest --quiet >"$log" 2>&1; then
  rm -f "$log"
  echo '{}'
  exit 0
fi

tail_output=$(tail -n 40 "$log")
rm -f "$log"
jq -n --arg out "$tail_output" '{
  followup_message: ("`./gradlew :tassel-core:jvmTest` failed after your changes. Fix the failures (do not edit testdata/golden or loosen tolerances), or explain why they are unrelated. Last output:\n\n" + $out)
}'
exit 0
