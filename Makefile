# Build and install Python package
SRC_HPP = $(wildcard src/*hpp)
SRC_CPP = $(wildcard src/*cpp)

src: $(SRC_HPP) $(SRC_CPP) meson.build
	@echo "****************************************************************"
	@echo "Build and install Python package"
	rm -rf build
	pip install .
	touch src
	@echo "Done"
	@echo "****************************************************************"

quick: $(SRC_HPP) $(SRC_CPP) meson.build
	@echo "****************************************************************"
	@echo "Build Python package locally"
	rm -rf build
	meson setup build
	meson compile -C build
	@echo "Done"
	@echo "****************************************************************"

# Generate Doxygen docs/html/*
docs: $(SRC_HPP) $(SRC_CPP) meson.build pyproject.toml README.md
	@echo "****************************************************************"
	@echo "Build Doxygen documentation"
	doxygen Doxyfile
	@echo "Done"
	@echo "****************************************************************"

# Sync dev branch with main
.PHONY: dev
dev:
	git checkout main
	git pull
	git checkout dev
	git merge main
	git push

# Clean-up
.PHONY: clean
clean:
	rm -rf build
	rm -rf docs/html