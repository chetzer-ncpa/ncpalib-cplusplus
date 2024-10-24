parts=math gtest
test-parts=$(addprefix test-, $(parts))
get_target = $(subst $(2)-,,$(1))
define test-lib =
$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS=$(INCLUDE_FLAGS)
endef

all: test

test: $(test-parts)

$(test-parts):
	$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS="$(INCLUDE_FLAGS)"  LIBRARY_FLAGS="$(LIBRARY_FLAGS)"


.PHONY: $(test-parts) $(parts)