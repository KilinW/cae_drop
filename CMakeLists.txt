cmake_minimum_required(VERSION 3.13)

set(PICO_BOARD pico_w)

include($ENV{PICO_SDK_PATH}/external/pico_sdk_import.cmake)

project(blink_w C CXX ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

pico_sdk_init()

add_executable(${PROJECT_NAME}
    main.cpp
)

# Add pico_stdlib library, add more if used
target_link_libraries(${PROJECT_NAME}
    pico_stdlib
    pico_cyw43_arch_none
    )

# enable usb output, disable uart output
pico_enable_stdio_usb(blink_w 1)
pico_enable_stdio_uart(blink_w 0)

# Need to generate UF2 file for upload to RP2040
pico_add_extra_outputs(${PROJECT_NAME})