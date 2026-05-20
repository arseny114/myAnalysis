# Минимальная версия CMake
cmake_minimum_required(VERSION 3.10)

# Название проекта
project(zh_invisible_analysis)

set(SRC_CPP ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Находим ROOT с нужными библиотеками
find_package(ROOT 6.30 CONFIG REQUIRED COMPONENTS 
    Core 
    Hist 
    Tree 
    RooFit 
    RooFitCore 
    RooStats 
    Physics
)

# Отладочная версия
set(CMAKE_BUILD_TYPE Debug)

# Путь к include/
include_directories(${PROJECT_SOURCE_DIR}/include)
include(${ROOT_USE_FILE})

# Создаем исполняемый файл
add_executable(${PROJECT_NAME} ${SRC_CPP})

# Линкуем рутовские библиотеки
target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        ROOT::Core 
        ROOT::Hist 
        ROOT::Tree 
        ROOT::RooFit 
        ROOT::RooFitCore 
        ROOT::RooStats 
        ROOT::Physics
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
