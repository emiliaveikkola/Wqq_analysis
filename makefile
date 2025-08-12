# ====== CONFIG ======
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXX      := g++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS) -Iinterface
LDFLAGS  := -shared $(ROOTLIBS)

SRC_DIR  := src
INC_DIR  := interface
OBJ_DIR  := obj
LIB_DIR  := lib
LIBNAME  := libAnalysis.so
LIBPATH  := $(LIB_DIR)/$(LIBNAME)

MACRO_DIR := make
ROOT     := root -l -b -q

# Pretty output (can be disabled with VERBOSE=1 to show full commands)
GREEN := \033[32m
RED   := \033[31m
BLUE  := \033[34m
BOLD  := \033[1m
DIM   := \033[2m
RESET := \033[0m

# Common ROOT exec prefix
ROOTLINE := $(ROOT) -e '(void)gSystem->Load("$(LIBPATH)")' -e 'gSystem->AddIncludePath("-Iinterface")' -e '.x '

# ====== DISCOVER SOURCES ======
SRCS := $(wildcard $(SRC_DIR)/*.cc)
OBJS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRCS))

# ====== JOBS (compile -> run) ======
JOBS := tagandprobe wqq strangejet

COMPILE_tagandprobe := $(MACRO_DIR)/mk_tpcompile.C
RUN_tagandprobe     := $(MACRO_DIR)/mk_TagandProbe.C

COMPILE_wqq := $(MACRO_DIR)/mk_Wqqcompile.C
RUN_wqq     := $(MACRO_DIR)/mk_Wqq.C

COMPILE_strangejet := $(MACRO_DIR)/mk_sjcompile.C
RUN_strangejet     := $(MACRO_DIR)/mk_StrangeJet.C

# ====== TOP-LEVEL TARGETS ======
.PHONY: all build list clean $(JOBS) run-% compile-%

all: $(JOBS)

# Build only the shared lib
build: $(LIBPATH)

# Each job: ensure lib exists, compile (if macro exists), then run
$(JOBS): %: build compile-% run-%

# ====== BUILD RULES ======
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIBPATH): $(OBJS)
	@mkdir -p $(LIB_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

# ====== MACRO RUNNERS ======
compile-%:
	@c="$(COMPILE_$*)"; \
	if [ -n "$$c" ] && [ -f "$$c" ]; then \
	  printf "$(BOLD)🔧  Compiling [%-12s]$(RESET) — %s\n" "$*" "$$c"; \
	  if [ "$(VERBOSE)" = "1" ]; then printf "$(DIM)%s %s$(RESET)\n" "$(ROOTLINE)" "$$c"; fi; \
	  $(ROOTLINE) "$$c" >/dev/null || { printf "$(RED)✗ compile failed$(RESET)\n"; exit 1; }; \
	  printf "$(GREEN)✓ compile ok$(RESET)\n"; \
	else \
	  printf "$(BLUE)ℹ︎  [%s] no compile macro found — skipping$(RESET)\n" "$*"; \
	fi

run-%:
	@r="$(RUN_$*)"; \
	if [ -z "$$r" ] || [ ! -f "$$r" ]; then \
	  printf "$(RED)!! [%s] run macro not found: %s$(RESET)\n" "$*" "$$r"; exit 1; \
	fi; \
	printf "$(BOLD)▶️   Running  [%-12s]$(RESET) — %s\n" "$*" "$$r"; \
	if [ "$(VERBOSE)" = "1" ]; then printf "$(DIM)%s %s$(RESET)\n" "$(ROOTLINE)" "$$r"; fi; \
	$(ROOTLINE) "$$r" || { printf "$(RED)✗ run failed$(RESET)\n"; exit 1; }; \
	printf "$(GREEN)✓ run ok$(RESET)\n"

list:
	@echo "Jobs:" $(JOBS)
	@for j in $(JOBS); do \
	  printf "  - %s\n      compile: %s\n      run:     %s\n" \
	    $$j "$$(eval echo \$$(echo COMPILE_$$j))" "$$(eval echo \$$(echo RUN_$$j))"; \
	done

clean:
	@echo "Cleaning…"
	@rm -rf $(OBJ_DIR) $(LIB_DIR)
	@find . -type f \( -name "*.so" -o -name "*.d" -o -name "*.pcm" \) -delete
	@echo "Done."