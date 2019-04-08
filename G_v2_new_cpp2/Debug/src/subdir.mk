################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Fix3DFunc.cpp \
../src/FixFunc.cpp \
../src/G_main_v2.cpp \
../src/Geval.cpp \
../src/vectoralgebra.cpp 

OBJS += \
./src/Fix3DFunc.o \
./src/FixFunc.o \
./src/G_main_v2.o \
./src/Geval.o \
./src/vectoralgebra.o 

CPP_DEPS += \
./src/Fix3DFunc.d \
./src/FixFunc.d \
./src/G_main_v2.d \
./src/Geval.d \
./src/vectoralgebra.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


