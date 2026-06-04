# Минимальная версия CMake
cmake_minimum_required(VERSION 3.10)

# Название проекта
project(zh_invisible_analysis)

set(SRC_CPP ${PROJECT_SOURCE_DIR}/src/main.cpp)
set(XGBOOST_ROOT "~/xgboost")

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
include_directories(${XGBOOST_ROOT}/include)
include(${ROOT_USE_FILE})

# Создаем исполняемый файл
add_executable(${PROJECT_NAME} ${SRC_CPP})

# Линкуем рутовские библиотеки и xgboost
target_link_directories(zh_invisible_analysis PRIVATE ${XGBOOST_ROOT}/lib)
target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        ROOT::Core 
        ROOT::Hist 
        ROOT::Tree 
        ROOT::RooFit 
        ROOT::RooFitCore 
        ROOT::RooStats 
        ROOT::Physics
        xgboost
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
