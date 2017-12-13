# コンパイラの指定　基本g++でいい
COMPILER  = g++
# OpenMPを使った並列化を行うときはここを有効にする
OPENMP  = -fopenmp
# WARNING = -wall
# DEBUG   = -pg
# コンパイルオプションの指定
CFLAGS    = -std=gnu++14 -O3 -funroll-loops -ffast-math -march=native $(OPENMP) $(WARNING)
# CFLAGS    = -std=gnu++14 -fopenmp
# リンクオプションの設定
ifeq "$(shell getconf LONG_BIT)" "64"
  LDFLAGS =
else
  LDFLAGS =
endif
# ライブラリを利用するときは入れる　（例:-lm）
LIBS      =
# ディレクトリの指定　それぞれヘッダーファイル/実行ファイル置き場/ソースコード置き場
INCLUDE   = -I./include
TARGET    = ./shock_tube
SRCDIR    = ./source
ifeq "$(strip $(SRCDIR))" ""
  SRCDIR  = .
endif

# SOURCES   = $(wildcard $(SRCDIR)/*.cpp) # 全部使ってコンパイルするとき
# コンパイルするファイルが決まっているときはshellコマンドで追加していく
SOURCES = $(shell ls $(SRCDIR)/shock_tube.cpp)\
$(shell ls $(SRCDIR)/riemann_function.cpp)\

# オブジェクトファイル置き場
OBJDIR    = ./object
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR  = .
endif
# ターゲットである実行ファイルが依存するオブジェクトファイルを指定
# SOURCESで指定したソースファイルから作成されたオブジェクトを全て使用することを想定
# substを使用して，SRC_DIRからOBJ_DIRに変換し，拡張子も.cppから.oにする
OBJECTS   = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
# 依存関係を示す中間ファイルの宣言
DEPENDS   = $(OBJECTS:.o=.d)

$(TARGET): $(OBJECTS) $(LIBS)
	$(COMPILER) $(CFLAGS)$(DEBUG) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

all: clean $(TARGET)

clean:
	-rm -f $(OBJECTS) $(DEPENDS) $(TARGET)
	-rm -f ./data/shock_tube_*.dat

# sourceファイルの依存関係が明記された.dファイルをインクルードする
# これによってheaderファイルのみが更新された場合でもmakeが実行されることになる
-include $(DEPENDS)
# ファイルを生成しないターゲットを明記
.PHONY: all clean
