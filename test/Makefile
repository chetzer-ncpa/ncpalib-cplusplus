export
parts=NCPA
test-parts=$(addprefix test-, $(parts))
get_target = $(subst $(2)-,,$(1))
# define test-lib =
# $(MAKE) -C $(call get_target,$@,test) test INCLUDE_FLAGS=$(INCLUDE_FLAGS)
# endef

# override env variables

#export BASE_DIR, INCLUDE_PATH_FLAGS, LIBRARY_PATH_FLAGS, STANDARD

all: test

test: $(test-parts)

$(test-parts):
	$(MAKE) -C $(call get_target,$@,test) clean
	$(MAKE) -C $(call get_target,$@,test) test
	$(MAKE) -C $(call get_target,$@,test) clean

.PHONY:  $(test-parts)