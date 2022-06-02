MKDIR = mkdir
CP = cp
RM = rm
MV = mv
TAR = tar
EXEC = dfvs

build:
	cmake -S . -B build/Release -DCMAKE_BUILD_TYPE=Release
	cmake --build build/Release
	$(MKDIR) -p dist
	$(CP) -f build/Release/$(EXEC) dist/
	@echo "Created: dist/$(EXEC)"

clean:
	@echo "Cleaning..."
	@$(RM) -rf build/*
	@$(RM) -rf dist/*
	@echo "Cleaning done."

publish:
	$(MKDIR) -p dist
	$(TAR) zcvf dist/$(EXEC).tgz  --exclude-from=.tarignore CMakeLists.txt src
	@echo "Created: dist/$(EXEC).tgz"

.PHONY: build clean publish
