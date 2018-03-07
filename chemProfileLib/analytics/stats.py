import math
# The implementation of these methods is originally from Mod-Seeker:
#
#Mod-seq data analysis pipeline
#Copyright (C) 2014  Yizhu Lin
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

def cmh(array):
    """ This function makes a
            Cochran-Mantel-Haenszel Chi-Squared Test for Count Data
            This test method is known in R, but unknown in scipy """

    reps = len(array) // 4
    array = [int(x) for x in array]
    chit1Sum = 0
    chit2Sum = 0
    ORt1Sum = 0
    ORt2Sum = 0
    SEt1_nSum = 0
    SEt1_dSum = 0
    SEt2_nSum = 0
    SEt2_dSum = 0
    SEt3_nSum = 0

    for i in range(0, reps):
        a = array[4 * i]
        b = array[4 * i + 1]
        c = array[4 * i + 2]
        d = array[4 * i + 3]
        n = float(a + b + c + d)

        n2 = n**2
        n3 = n**3
        n3m2 = n3 - n2
        chit1 = a - (a + b) * (a + c) / n
        chit2 = (a + b) * (a + c) * (b + d) * (c + d) / n3m2

        chit1Sum += chit1
        chit2Sum += chit2

        ORt1 = a * d / n
        ORt2 = b * c / n
        ORt1Sum += ORt1
        ORt2Sum += ORt2

        SEt1_n = (a + d) * a * d / n2
        SEt1_d = a * d / n
        SEt2_n = (b + c) * b * c / n2
        SEt2_d = b * c / n
        SEt3_n = ((a + d) * b * c + (b + c) * a * d) / n2
        SEt1_nSum += SEt1_n
        SEt1_dSum += SEt1_d
        SEt2_nSum += SEt2_n
        SEt2_dSum += SEt2_d
        SEt3_nSum += SEt3_n

    chi = (abs(chit1Sum) - 0.5) ** 2 / chit2Sum
    OR = ORt1Sum / ORt2Sum
    var = SEt1_nSum / (2 * SEt1_dSum ** 2) + SEt2_nSum / (2 * SEt2_dSum ** 2) + SEt3_nSum / (
        2 * SEt1_dSum * SEt2_dSum)
    SE = math.sqrt(var)
    z = 1.959964
    ORL = OR * math.exp(SE * (-z))
    ORU = OR * math.exp(SE * z)

    pvalue = pchisq(chi)

    return (chi, pvalue, OR, ORL, ORU)

def pchisq(x):
    #return p-value of chi-sq test with df=1
    return 1 - math.erf(math.sqrt(0.5*float(x)))

def BHcontrol(pvalues, FDR = 0.05):
    pvalues.sort()
    m = float(len(pvalues))
    p = 0
    i = 0

    while p <= ((i + 1) / m) * FDR and i < m:
        p = pvalues[i]
        i += 1

    p_sig = pvalues[max(i - 1, 0)]

    if p_sig <= FDR:
        # print "BHcontrol:", p_sig, m, i-1
        return p_sig
    else:
        return None