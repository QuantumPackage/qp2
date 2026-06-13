default: config/Makefile
	bash -c "source quantum_package.rc ; $(MAKE) -f config/Makefile"

config/Makefile:
	@bash -c ' echo '' ; echo xxxxxxxxxxxxxxxxxx ; echo "QP is not configured yet. Please run the ./configure command" ; echo xxxxxxxxxxxxxxxxxx ; echo ''  ; ./configure --help' | more
