# -*- coding: UTF8 -*-

#############################################################################
# Author: Guillaume Bouvier -- guillaume.bouvier@pasteur.fr                 #
# https://research.pasteur.fr/en/member/guillaume-bouvier/                  #
# Copyright (c) 2022 Institut Pasteur                                       #
#                                                                           #
#                                                                           #
#  Redistribution and use in source and binary forms, with or without       #
#  modification, are permitted provided that the following conditions       #
#  are met:                                                                 #
#                                                                           #
#  1. Redistributions of source code must retain the above copyright        #
#  notice, this list of conditions and the following disclaimer.            #
#  2. Redistributions in binary form must reproduce the above copyright     #
#  notice, this list of conditions and the following disclaimer in the      #
#  documentation and/or other materials provided with the distribution.     #
#  3. Neither the name of the copyright holder nor the names of its         #
#  contributors may be used to endorse or promote products derived from     #
#  this software without specific prior written permission.                 #
#                                                                           #
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS      #
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT        #
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR    #
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT     #
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT         #
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,    #
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY    #
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      #
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE    #
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.     #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#                                                                           #
#############################################################################

import sys
from datetime import datetime

class Progress:
    def __init__(self, n_step, delta = 1, label = None):
        """

        • n_step: total number of step

        • delta: delta in percent (default 1%)

        • label: Optional string to display

        """
        self.n_step = n_step
        self.progress = set([ int(self.n_step*p/100) for p in range(0,100,delta)[1:] ])
        self.progress.update([self.n_step,])
        self.c = 0
        self.delta = delta
        self.t1 = datetime.now()
        self.label = label

    def count(self, report=None):
        """
        • report: A string or value you want to report
        """
        self.c += 1
        if self.c in self.progress or self.c == 1:
            t2 = datetime.now()
            percent = float(self.c)*100/(self.n_step)
            if self.c > 1:
                eta = (t2 - self.t1) * int((100 - percent) / self.delta)
            else:
                eta = (t2 - self.t1) * int((100 - percent) / percent)
            if self.label is None:
                string = "%d %% ETA: %s"%(percent, eta)
            else:
                string = "%s: %d %% ETA: %s"%(self.label, percent, eta)
            if report is not None:
                string+=" | %s"%report
            print(string)
            sys.stdout.flush()
            self.t1 = t2
