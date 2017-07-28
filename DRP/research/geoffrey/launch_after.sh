#!/bin/bash
while [ -d /proc/$1 ] ; do
    sleep 100
done && $2