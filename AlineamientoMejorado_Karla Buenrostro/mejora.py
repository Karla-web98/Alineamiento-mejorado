def align_sequences(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_open=-3, gap_extend=-1):
    n = len(seq1)
    m = len(seq2)

    # Inicialización de matrices
    score = [[0] * (m + 1) for _ in range(n + 1)]
    gap1 = [[0] * (m + 1) for _ in range(n + 1)]  # Gaps verticales (en seq1)
    gap2 = [[0] * (m + 1) for _ in range(n + 1)]  # Gaps horizontales (en seq2)

    for i in range(1, n + 1):
        score[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        score[0][j] = gap_open + (j - 1) * gap_extend

    # Llenado de la matriz
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)

            gap_seq1 = max(
                score[i - 1][j] + gap_open,
                gap1[i - 1][j] + gap_extend
            )
            gap_seq2 = max(
                score[i][j - 1] + gap_open,
                gap2[i][j - 1] + gap_extend
            )

            score[i][j] = max(match, gap_seq1, gap_seq2)
            gap1[i][j] = gap_seq1
            gap2[i][j] = gap_seq2

    return score[n][m]


# =======================
# PRUEBAS
# =======================

if __name__ == "__main__":
    ejemplos = [
        ("AGCT", "AGT"),
        ("GATTACA", "GCATGCU"),
        ("ACTGAC", "ACTG---GAC")
    ]

    for s1, s2 in ejemplos:
        fitness = align_sequences(s1.replace("-", ""), s2.replace("-", ""))
        print(f"Secuencia 1: {s1}")
        print(f"Secuencia 2: {s2}")
        print(f"Fitness obtenido: {fitness}")
        print("-" * 40)
