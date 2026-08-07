#!/bin/sh
until [ ! -z "$(grep -c 'convertible' /tmp/claude-0/-home-user-integer-solver/7ccf248f-5604-592b-968f-0791d039bf11/tasks/blx1ason5.output 2>/dev/null)" ] && [ "$(grep -c 'convertible' /tmp/claude-0/-home-user-integer-solver/7ccf248f-5604-592b-968f-0791d039bf11/tasks/blx1ason5.output)" -ge 2 ]; do sleep 5; done
