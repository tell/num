
all:
	$(MAKE) -C src
	$(MAKE) -C test
	$(MAKE) -C bench

clean:
	$(MAKE) -C bench clean
	$(MAKE) -C test clean
	$(MAKE) -C src clean
