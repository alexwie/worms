#!/bin/sh

H="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0"

(for h in $H; do ./worms -L 4 -T 0.25 -H $h; done) | awk '$1=="Magnetic" { h=$4 } $1=="Uniform" && $2=="Magnetization" {print h, $4, $6}'
