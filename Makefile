default: build.ninja
	bash -c "source quantum_package.rc ; ninja"

build.ninja:
	@bash -c ' echo '' ; echo xxxxxxxxxxxxxxxxxx ; echo "QP is not configured yet. Please run the ./configure command" ; echo xxxxxxxxxxxxxxxxxx ; echo ''  ; ./configure --help' | more
