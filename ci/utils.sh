#!/bin/bash

# Various utility functions used through CI.

# Finds Cargo's `OUT_DIR` directory from the most recent build.
#
# This requires one parameter corresponding to the target directory
# to search for the build output.
cargo_out_dir() {
    # This works by finding the most recent stamp file, which is produced by
    # every vartrix build.
    target_dir="$1"
    find "$target_dir" -name vartrix-stamp -print0 \
      | xargs -0 ls -t \
      | head -n1 \
      | xargs dirname
}

host() {
    case "$TRAVIS_OS_NAME" in
        linux)
            case "${TRAVIS_CPU_ARCH:-amd64}" in
                arm64)
                    echo aarch64-unknown-linux-gnu
                    ;;
                *)
                    echo x86_64-unknown-linux-gnu
                    ;;
            esac
            ;;
        osx)
            case "${TRAVIS_CPU_ARCH:-amd64}" in
                arm64)
                    echo aarch64-apple-darwin
                    ;;
                *)
                    echo x86_64-apple-darwin
                    ;;
            esac
            ;;
    esac
}

architecture() {
    case "$TARGET" in
        x86_64-*)
            echo amd64
            ;;
        i686-*|i586-*|i386-*)
            echo i386
            ;;
        aarch64-*)
            echo aarch64
            ;;
        arm*-unknown-linux-gnueabihf)
            echo armhf
            ;;
        *)
            die "architecture: unexpected target $TARGET"
            ;;
    esac
}

gcc_prefix() {
    case "$(architecture)" in
        armhf)
            echo arm-linux-gnueabihf-
            ;;
        *)
            return
            ;;
    esac
}

is_ssse3_target() {
    case "$(architecture)" in
        amd64) return 0 ;;
        *)     return 1 ;;
    esac
}

is_x86() {
    case "$(architecture)" in
      amd64|i386) return 0 ;;
      *)          return 1 ;;
    esac
}

is_arm() {
    case "$(architecture)" in
        armhf) return 0 ;;
        *)     return 1 ;;
    esac
}

is_aarch64() {
    case "$(architecture)" in
        aarch64) return 0 ;;
        *)       return 1 ;;
    esac
}

is_linux() {
    case "$TRAVIS_OS_NAME" in
        linux) return 0 ;;
        *)     return 1 ;;
    esac
}

is_osx() {
    case "$TRAVIS_OS_NAME" in
        osx) return 0 ;;
        *)   return 1 ;;
    esac
}
