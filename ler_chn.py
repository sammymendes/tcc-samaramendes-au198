#!/usr/bin/env python3
"""
Leitor de espectros ORTEC no formato binario .Chn (integer format).

Uso:
    python3 ler_chn.py arquivo.Chn [saida.csv]

Layout do arquivo (little-endian):
  Cabecalho (32 bytes)
    0   int16   -1 (assinatura do formato)
    2   int16   numero do detector
    4   int16   numero do segmento
    6   char[2] segundos do inicio da contagem
    8   int32   tempo real   (unidades de 20 ms)
    12  int32   tempo vivo   (unidades de 20 ms)
    16  char[8] data de inicio DDMMMYY + flag de seculo ('1' => 20xx)
    24  char[4] hora de inicio HHMM
    28  int16   canal inicial
    30  int16   numero de canais
  Espectro: N x int32 (contagens por canal)
  Rodape (512 bytes, opcional)
    0   int16   -101 ou -102
    4   float32 calibracao em energia: termo constante, ganho, termo quadratico
    16  float32 calibracao de forma de pico (3 termos)
    28  char[4] unidade de energia (ex.: 'keV')
    256 pstring descricao do detector  (1 byte de tamanho + 63 bytes)
    320 pstring descricao da amostra   (1 byte de tamanho + 63 bytes)
"""

import struct
import sys

MESES = {m: i + 1 for i, m in enumerate(
    ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
     "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"])}


def _pstring(buf, offset):
    n = buf[offset]
    return buf[offset + 1: offset + 1 + n].decode("ascii", "replace").strip()


def ler_chn(caminho):
    with open(caminho, "rb") as fh:
        d = fh.read()

    assinatura, detector, segmento = struct.unpack_from("<hhh", d, 0)
    if assinatura != -1:
        raise ValueError("Assinatura invalida: nao parece um .Chn ORTEC")

    segundos = d[6:8].decode("ascii", "replace")
    tempo_real, tempo_vivo = struct.unpack_from("<ii", d, 8)
    data_txt = d[16:24].decode("ascii", "replace")
    hora_txt = d[24:28].decode("ascii", "replace")
    canal_inicial, n_canais = struct.unpack_from("<hh", d, 28)

    contagens = list(struct.unpack_from("<%di" % n_canais, d, 32))

    dia = int(data_txt[0:2])
    mes = MESES.get(data_txt[2:5].upper(), 0)
    ano = int(data_txt[5:7]) + (2000 if data_txt[7:8] == "1" else 1900)
    inicio = "%02d/%02d/%04d %s:%s:%s" % (
        dia, mes, ano, hora_txt[0:2], hora_txt[2:4], segundos)

    info = {
        "arquivo": caminho,
        "detector": detector,
        "segmento": segmento,
        "inicio_contagem": inicio,
        "tempo_real_s": tempo_real * 0.02,
        "tempo_vivo_s": tempo_vivo * 0.02,
        "canal_inicial": canal_inicial,
        "n_canais": n_canais,
        "contagens": contagens,
        "calib_energia": None,
        "calib_forma": None,
        "unidade_energia": None,
        "descricao_detector": None,
        "descricao_amostra": None,
    }

    tempo_morto = tempo_real - tempo_vivo
    info["tempo_morto_pct"] = 100.0 * tempo_morto / tempo_real if tempo_real else 0.0

    rodape = 32 + 4 * n_canais
    if len(d) >= rodape + 32:
        f = d[rodape:]
        tipo = struct.unpack_from("<h", f, 0)[0]
        if tipo in (-101, -102):
            info["calib_energia"] = struct.unpack_from("<3f", f, 4)
            info["calib_forma"] = struct.unpack_from("<3f", f, 16)
            info["unidade_energia"] = f[28:32].split(b"\x00")[0].decode("ascii", "replace")
            if len(f) >= 384:
                info["descricao_detector"] = _pstring(f, 256)
                info["descricao_amostra"] = _pstring(f, 320)

    if info["calib_energia"]:
        a, b, c = info["calib_energia"]
        info["energias"] = [a + b * ch + c * ch * ch
                            for ch in range(canal_inicial, canal_inicial + n_canais)]
    else:
        info["energias"] = None

    return info


def para_csv(info, saida):
    with open(saida, "w", encoding="utf-8") as fh:
        if info["energias"]:
            fh.write("canal,energia_%s,contagens\n" % (info["unidade_energia"] or "keV"))
            for i, (e, c) in enumerate(zip(info["energias"], info["contagens"])):
                fh.write("%d,%.4f,%d\n" % (info["canal_inicial"] + i, e, c))
        else:
            fh.write("canal,contagens\n")
            for i, c in enumerate(info["contagens"]):
                fh.write("%d,%d\n" % (info["canal_inicial"] + i, c))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    info = ler_chn(sys.argv[1])
    print("Amostra ............ %s" % info["descricao_amostra"])
    print("Detector ........... %s (n. %d)" % (info["descricao_detector"], info["detector"]))
    print("Inicio da contagem . %s" % info["inicio_contagem"])
    print("Tempo real / vivo .. %.2f s / %.2f s (tempo morto %.2f%%)" % (
        info["tempo_real_s"], info["tempo_vivo_s"], info["tempo_morto_pct"]))
    print("Canais ............. %d (a partir do canal %d)" % (info["n_canais"], info["canal_inicial"]))
    if info["calib_energia"]:
        a, b, c = info["calib_energia"]
        print("Calibracao ......... E = %.6g + %.6g*ch + %.6g*ch^2  [%s]" % (
            a, b, c, info["unidade_energia"]))
    print("Contagens totais ... %d" % sum(info["contagens"]))
    if len(sys.argv) > 2:
        para_csv(info, sys.argv[2])
        print("CSV gravado em ..... %s" % sys.argv[2])
