parts=NCPA
test-parts=$(addprefix test-, $(parts))
get_target = $(subst $(2)-,,$(1))
define test-lib =
#$(MAKE) -C $(call get_target,$@,test) clean
$(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS=$(INCLUDE_FLAGS)
endef

all: test

test: $(test-parts)

$(test-parts):
	$(test-lib)

.PHONY:  $(test-parts)