################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../alglib/alglibinternal.cpp \
../alglib/alglibmisc.cpp \
../alglib/ap.cpp \
../alglib/dataanalysis.cpp \
../alglib/diffequations.cpp \
../alglib/fasttransforms.cpp \
../alglib/integration.cpp \
../alglib/interpolation.cpp \
../alglib/linalg.cpp \
../alglib/optimization.cpp \
../alglib/solvers.cpp \
../alglib/specialfunctions.cpp \
../alglib/statistics.cpp 

OBJS += \
./alglib/alglibinternal.o \
./alglib/alglibmisc.o \
./alglib/ap.o \
./alglib/dataanalysis.o \
./alglib/diffequations.o \
./alglib/fasttransforms.o \
./alglib/integration.o \
./alglib/interpolation.o \
./alglib/linalg.o \
./alglib/optimization.o \
./alglib/solvers.o \
./alglib/specialfunctions.o \
./alglib/statistics.o 

CPP_DEPS += \
./alglib/alglibinternal.d \
./alglib/alglibmisc.d \
./alglib/ap.d \
./alglib/dataanalysis.d \
./alglib/diffequations.d \
./alglib/fasttransforms.d \
./alglib/integration.d \
./alglib/interpolation.d \
./alglib/linalg.d \
./alglib/optimization.d \
./alglib/solvers.d \
./alglib/specialfunctions.d \
./alglib/statistics.d 


# Each subdirectory must supply rules for building sources it contributes
alglib/%.o: ../alglib/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


