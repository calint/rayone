################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/rayone.cpp 

OBJS += \
./src/rayone.o 

CPP_DEPS += \
./src/rayone.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	clang++ -I/System/Library/Frameworks/GLUT.framework/Versions/Current/Headers/ -I/System/Library/Frameworks/OpenGL.framework/Versions/Current/Headers/ -O0 -g3 -pedantic -pedantic-errors -Wall -Wextra -Werror -Wconversion -c -fmessage-length=0 -std=gnu++11 -Wfatal-errors -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


