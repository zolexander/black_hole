#!/bin/bash

# --- Farbdefinitionen für Ausgaben --- #
RED='\033[;31m'
YELLOW='\033[1;33m'
NC='\033[0m'
MAGENTA='\033[0;35m'

SCRIPT_NAME=$0
CMAKELISTS="CMakeLists.txt"

# --- Hilfefunktion --- #
display_help() {
    echo -e "Usage: ${SCRIPT_NAME} [-b|--build DIRECTORY] [-c] [-h|--help]\n"
    echo "Options:"
    echo -e "  ${YELLOW}-b, --build DIRECTORY      ${NC}Specify a directory for building"
    echo -e "  ${YELLOW}-c, --clear-build          ${NC}Clear the complete build directory content with cache together"
    echo -e "  ${YELLOW}-h, --help                 ${NC}Display this help message"
    exit 0
}

# --- Option Parsing --- #
OPTS=$(getopt -o b:ch --long build:,clear-build,help -n "${SCRIPT_NAME}" -- "$@")
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to parse options${NC}\n" >&2
    display_help
    exit 1
fi
eval set -- "$OPTS"

# --- Hilfsfunktionen --- #
make_colored() {
    make 2>&1 | sed \
        -e "s/\(error:\)/$(printf '\033[1;31m')\1$(printf '\033[0m')/g" \
        -e "s/\(warning:\)/$(printf '\033[1;33m')\1$(printf '\033[0m')/g"
}

parse_cmake_executables() {
    local result=()
    while IFS= read -r line; do
        if [[ $line == *"add_executable("* ]]; then
            local exec=${line#*(}
            if [[ $exec == *".cpp"* ]]; then
                result+=("${exec% *}")
            else
                result+=("${exec}")
            fi
        fi
    done <"../${CMAKELISTS}"
    echo "${result[@]}"
}

delete_executables() {
    for exe in "$@"; do
        echo -e "${RED}./build/$exe${NC}"
        rm -f "$exe"
    done
}

# --- Standardwerte --- #
BUILD_DIRECTORY=""
CLEAR_BUILD=0

# --- Optionen auswerten --- #
while true; do
    case "$1" in
        -b|--build)
            BUILD_DIRECTORY="$2"
            shift 2
            ;;
        -c|--clear-build)
            CLEAR_BUILD=1
            shift
            ;;
        -h|--help)
            display_help
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            display_help
            exit 1
            ;;
    esac
done

# --- Vorbedingungen prüfen --- #
if [ -z "$BUILD_DIRECTORY" ]; then
    echo -e "${RED}No output directory was given.${NC}"
    exit 1
fi
if [ -f "$BUILD_DIRECTORY" ]; then
    echo -e "${RED}Error: output directory is not a directory.${NC}"
    exit 1
fi
if [ ! -f "$CMAKELISTS" ]; then
    echo -e "${RED}The CMakeLists.txt file was not found in the project directory. Exiting...${NC}"
    exit 1
fi
if [ -f "CMakeCache.txt" ]; then
    echo "There is a CMakeCache.txt file in the project root directory. I will delete it before compiling"
    rm -rf "CMakeCache.txt"
fi
if [ ! -d "$BUILD_DIRECTORY" ]; then
    mkdir -p "$BUILD_DIRECTORY"
fi

# --- In das Build-Verzeichnis wechseln --- #
cd "$BUILD_DIRECTORY" || exit 1

# --- Build-Verzeichnis bereinigen --- #
if [ "$(ls -A)" ] && [ "$CLEAR_BUILD" -eq 1 ]; then
    echo -e "${YELLOW}Build directory is not empty and the CLEAR_BUILD flag is set. Clearing...${NC}"
    rm -rf *
else
    echo "The $BUILD_DIRECTORY is not empty. I will delete only the executable files."
    rm -rf CMakeCache.txt
    rm -rf MakeFile
    executables=($(parse_cmake_executables))
    if [ "${#executables[@]}" -gt 0 ]; then
        echo "The following files will be deleted:"
        delete_executables "${executables[@]}"
    fi
fi

# --- Build-Prozess --- #
echo -e "${YELLOW}"
cmake ..
echo -e "${NC}"
if [ $? -eq 0 ]; then
    echo -e "Making project...\n${MAGENTA}"
    make_colored -j "$(nproc)"
    if [ $? -eq 0 ]; then
        echo -e "${NC}${YELLOW}The file has been successfully compiled.${NC}"
        exit 0
    else
        echo -e "${NC}${RED}Some error occurred during make.${NC}"
        exit 1
    fi
else
    echo -e "${RED}Some error occurred during the CMake process.${NC}"
    exit 1
fi