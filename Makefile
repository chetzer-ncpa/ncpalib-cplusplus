parts=test
test-parts=$(addprefix test-, $(parts))
get_target = $(subst $(2)-,,$(1))
BASE_DIR:=$(strip $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST)))))
INCLUDE_FLAGS=-I$(BASE_DIR)/include
define test-target =
#$(MAKE) -C $(call get_target,$@,test) clean
$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS=$(INCLUDE_FLAGS)
endef

all: test

test: $(test-parts)

$(test-parts):
	$(test-target)

.PHONY:  $(test-parts)