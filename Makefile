.PHONY: all
.PHONY: clean
.PHONY: test

all: redsvd

clean:
	rm -r bin include lib redsvd_build

redsvd:
	./waf configure --blddir redsvd_build --prefix .
	./waf
	./waf install

test: redsvd
	bin/redsvd -i test_files/matrix.txt -o output -f sparse
