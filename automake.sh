#!/bin/bash

RED='\033[;31m'
YELLOW='\033[1;33m'
NC='\033[0m'
MAGENTA='\033[0;35m'
SCRIPT_NAME=$0
displayhelp() {

    echo -e "Usage: ${SCRIPT_NAME} [-b|--build DIRECTORY] [-c] [-h|--help]\n"
    echo "Options:"
    echo -e "  ${YELLOW}-b, --build DIRECTORY      ${NC}Specify a directory for building"
    echo -e "  ${YELLOW}-c, --clear-build          ${NC}Clear the complete build directory content with cache together"
    echo -e "  ${YELLOW}-h, --help                 ${NC}Display this help message"
    exit 0
}
OPTS=$(getopt -o b:ch --long build:,directory:,clear-build,help -n "${SCRIPT_NAME}" -- "$@")
if [ $? -ne 0 ]; then
	echo -e "${RED}Failed to parse options${NC}\n" >&2
	displayhelp
	exit 1
fi


eval set -- "$OPTS"
make_colored() {
	make 2>&1 | sed \
    		-e "s/\(error:\)/$(printf '\033[1;31m')\1$(printf '\033[0m')/g" \
    		-e "s/\(warning:\)/$(printf '\033[1;33m')\1$(printf '\033[0m')/g"


}
BUILD_DIRECTORY=""
CLEARBUILD=0
CMAKELISTS="CMakeLists.txt"
while true; do
	case "$1" in
		-b | --build)
		BUILD_DIRECTORY="$2"
		shift 2
		;;
	 	-h | --help)
		displayhelp
		shift
		;;
		-c | --clear-build)
		CLEARBUILD=1
		shift
		;;
		--)
		shift
		break
		;;
		*)
		displayhelp
		exit 1
		;;
	esac
done
## Display the results
if [ "$BUILD_DIRECTORY" =  '' ]; then
	echo -e "${RED}No output directory was given.${NC}"
	exit 0
fi
if [ -f "$BUILD_DIRECTORY" ]; then
	echo -e "${RED}Error: output directory is not a directory.${NC}"
	exit 1
fi

if [ ! -f "$CMAKELISTS" ]; then
	echo -e "${RED}The CMakeLists.txt file was not found in the project directory. Exiting..."
	exit 1
fi
if [ -f "CMakeCache.txt" ]; then
	echo "There is a CMakeCache.txt file in the project root directory. I will delete it before compiling"
	rm -rf "CMakeCache.txt"
fi
if [ ! -d "$BUILD_DIRECTORY" ]; then
	mkdir -p "$BUILD_DIRECTORY"
fi
declare -a RESULT=()

parsecmake() {
	while IFS= read -r line; do
	if [[ $line == *"add_executable("* ]]; then
		EXECUTABLE=${line#*(}
		if [[ $EXECUTABLE == *".cpp"* ]]; then
		RESULT+=("${EXECUTABLE% *}")
		else
		RESULT+=("${EXECUTABLE}")
		fi
	fi
	done <"../${CMAKELISTS}"
	echo "${RESULT[@]}"
}
deleteFiles() {
   for i in "${RESULT[@]}";
      do
	  	  echo -e "${RED}./build/$i${NC}"
          rm -f "$i"
      done
}
cd "$BUILD_DIRECTORY" || return

if [ "$(ls -A)" ] && [ "${CLEARBUILD}" -eq 1 ];then
	echo -e "${YELLOW}Build directory is not empty and the CLEARBUILD flag is setted to true.\n I clear the cache before compiling...${NC}"
	rm -rf *
else
	echo "The $BUILD_DIRECTORY is not empty. I will delete only the executable files."
	rm -rf CMakeCache.txt
	rm -rf MakeFile
	parsecmake
	echo "The following files will be deleted:"
	deleteFiles
fi
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
		echo  -e"${NC}${RED}Some error concurred.${NC}"
			exit 1
	fi
else
	echo -e "${RED}Some error was encourred at the Cmake process.\n"
fi

