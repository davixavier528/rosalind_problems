def countnucleotides(DNA):
    countA = DNA.count("A")
    countC = DNA.count("C")
    countG = DNA.count("G")
    countT = DNA.count("T")
    print(countA, countC, countG, countT)

countnucleotides('AATGACTCTTACCTGCCGGCAGGTTTTAGACTCAATAAGTGATTGGAAAACTGGACCTTGATGGATTCGGTAGGGAACCAATGCTAAGTCACGTTAAGGACAAGTGACATGACTTTGTACTCTTCGTAACGAGCCTGTTGCACGCGTAGAGTGGTGTCGCTCAGTCCGGCAGACGGCGTGGCATAAATTGTAACCATCAAGGTTTTAGCTCAGAACATTGCTTACATCACGTTTACTCTTAACTGGCTAAATACCTTGTATGATAGACGCCCCAGCAAAGCCACAGAGAAGACCAACGCGGCACCAACCAATGATGTAGGGATCATGAGAGCCGGTTGATCAGGAAATAGCTAGGAGTGTGTAACGTTATTGCAGGACCGGGTACATGTACCGACGATGTCCTACGTGAAAATCGCTGCCTTCACACCTAGTGCAGGATAGCTGTCGCTCGATGGATTGATGGTTTGAGCGGCCTGCTATGATTCAGTACTCAGCAGGCACTTTCAGGCACAAAAAGGCAGCGGGAAATAGGGAAATCCCATGTTATCAGTGTTGCCTCCTAGCCCAACAGTCTAAAGGACAGATGAGGACACCCGCACTGCATAAGGCGCGGGCCGCGCAACGGTCCACTCATTCATCAAAGCGAAGCCCTGCATTGCTCGAGCGTTACGGCGGTCCTGAATTATTGAGGGAATGTAGGCCAGCGTCGTCGGCACACAGTGCCGAGCCGGGTAGGCCGGGAAACGCAGTTTGGTTCGCTGAGCTGCATACACGCGTACTCGCGGCGTAGACACCAGCGAGGTTTGTTGATATTAATACCTCCCCCTACGCTATCCAGATCATTGTTCCGGCATCCGAATTCGGCTCTGTCCACGTACACACCCTGTAGGGATCCTTGCGGTTCCCTCGATTAGGTTGGAAAGTAGTCAGGCTACCGCTGCTTTAGCAGACTTTTG')