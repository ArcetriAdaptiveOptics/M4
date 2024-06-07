#this module contains the configuration parameters for collecting and analysizing the IFF

#notes:
# the current shape is intended to allow a first startup of the procedure.
# make sure to re-shape the module into a more "formal" scheme (e.g. class, yaml, ...)

#### --- convert this file into a yaml and provide libraries to read it. so you can copy into the data folder and re-read it to analyze the dataset

# the following part contains the parameters for the initial padding and triggering sequence
numberOfTriggerZeros = 25 #number of empty (0 ampl.) frames before the trigger command

triggerModeId = 6 #index of the mirror mode to be used as a template for the triggering sequence
triggerModeAmp = 5e-6 #amplitude of the trigger command


#the following section contains the parameters for the registration sequence
numberOfRegistrationZeros = 1 #number of empty (0 ampl.) frames at the beginning of the registration sequence
registrationAmp = 500e-9 #amplitude of the registration command
registrationTemplate = (1,-1)
registrationCommandId = (110, 290,585,609,745) #actuator Ids to be used as references for the registration, in segment-order


