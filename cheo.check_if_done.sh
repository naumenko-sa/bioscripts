#!/bin/bash

if grep -q "END:" $1;
then
    echo $1 "done"
fi 