parts=test
test-parts=$(addprefix test-, $(parts))
get_target = $(subst $(2)-,,$(1))
BASE_DIR:=$(strip $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST)))))
INCLUDE_FLAGS:=-I$(BASE_DIR)/include -I/Users/claus/miniconda3/include
LIBRARY_FLAGS:=-L/Users/claus/miniconda3/lib -rpath /Users/claus/miniconda3/lib
define test-target =
$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS="$(INCLUDE_FLAGS)"
endef

all: test

test: $(test-parts)

$(test-parts):
	$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS="$(INCLUDE_FLAGS)" LIBRARY_FLAGS="$(LIBRARY_FLAGS)"

.PHONY:  test $(test-parts)
