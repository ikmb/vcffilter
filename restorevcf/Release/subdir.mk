################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../RestoreArgs.cpp \
../restorevcf.cpp 

CPP_DEPS += \
./RestoreArgs.d \
./restorevcf.d 

OBJS += \
./RestoreArgs.o \
./restorevcf.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean--2e-

clean--2e-:
	-$(RM) ./RestoreArgs.d ./RestoreArgs.o ./restorevcf.d ./restorevcf.o

.PHONY: clean--2e-

