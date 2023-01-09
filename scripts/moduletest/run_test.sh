#!/bin/bash

# FACTS main Directory
cd ../

# Set log file paths
LOG_STDOUT=$TESTSCRIPT_DIR/test/${MODULE_SET}.${MODULE}/test.out
LOG_STDERR=$TESTSCRIPT_DIR/test/${MODULE_SET}.${MODULE}/test.err

echo "Launching test script.
Logging output to :
$LOG_STDOUT
$LOG_STDERR."
