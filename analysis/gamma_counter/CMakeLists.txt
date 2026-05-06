# Минимальная версия CMake
cmake_minimum_required(VERSION 2.8)

# Название проекта
project(gamma_counter)

set(SRC_CPP ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Находим ROOT
find_package(ROOT 6.30 CONFIG REQUIRED)

# Отладочная версия
set(CMAKE_BUILD_TYPE Debug)

# Путь к include/
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# Создаем исполняемый файл
add_executable(${PROJECT_NAME} ${SRC_CPP})

# Линкуем рутовские библиотеки
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})
