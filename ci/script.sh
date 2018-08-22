#!/bin/bash

# build, test and generate docs in this phase

set -ex

. "$(dirname $0)/utils.sh"

main() {
    # Test a normal debug build.
    cargo build --target "$TARGET" --verbose --all

    # Show the output of the most recent build.rs stderr.
    set +x
    stderr="$(find "target/$TARGET/debug" -name stderr -print0 | xargs -0 ls -t | head -n1)"
    if [ -s "$stderr" ]; then
      echo "===== $stderr ====="
      cat "$stderr"
      echo "====="
    fi
    set -x

    # sanity check the file type
    file target/"$TARGET"/debug/vartrix

    # Run tests
    cargo test --target "$TARGET" --verbose

}

main
