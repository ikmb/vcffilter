################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../RemoveArgs.cpp \
../removesamples.cpp 

CPP_DEPS += \
./RemoveArgs.d \
./removesamples.d 

OBJS += \
./RemoveArgs.o \
./removesamples.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean--2e-

clean--2e-:
	-$(RM) ./RemoveArgs.d ./RemoveArgs.o ./removesamples.d ./removesamples.o

.PHONY: clean--2e-

