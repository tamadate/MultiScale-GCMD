NAME        = MD.out
SRCS        = ./src/*cpp ./src/*/*cpp ./src/*/*/*cpp
OBJS        = $(SRCS:.cpp=.o)
CXX         = mpicxx
CXXFLAGS    = -std=c++11 -fopenmp

.PHONY: all
all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) -O3 $(CXXFLAGS) -o $(NAME) $(OBJS)

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: fclean
fclean: clean
	$(RM) $(NAME)

.PHONY: re
re: fclean all
