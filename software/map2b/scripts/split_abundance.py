#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MAP2B-Cross-domain 丰度表拆分脚本
按 Kingdom 把 Abundance.xls 拆成两张表并各自重新归一化（每个样品列加和=1）：
  - 微生物表 Abundance.micro.xls  : Kingdom in (Bacteria, Archaea, Fungi)
  - 宏观生物表 Abundance.macro.xls: 其余 Kingdom（动物等）
默认删除总表 Abundance.xls；加 --keep-total 保留。
用法: python3 split_abundance.py -i <Abundance.xls> -o <输出目录> [--keep-total]
"""
import argparse, os, sys

MICRO_KINGDOMS = {'Bacteria', 'Archaea', 'Fungi'}

def read_abundance(path):
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        rows = []
        for line in f:
            t = line.rstrip('\n').split('\t')
            if len(t) >= len(header):
                rows.append(t)
    return header, rows

def write_table(path, header, rows, samples):
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')
    # 校验每样品列加和=1
    for si in samples:
        s = sum(float(r[si]) for r in rows if r[si].strip())
        print(f'  [check] {path} 列 {header[si]}: 加和={round(s, 6)}')

def main():
    parser = argparse.ArgumentParser(description='Split Abundance.xls into micro/macro tables')
    parser.add_argument('-i', dest='input', required=True)
    parser.add_argument('-o', dest='outdir', required=True)
    parser.add_argument('--keep-total', dest='keep_total', action='store_true')
    args = parser.parse_args()

    header, rows = read_abundance(args.input)
    samples = list(range(7, len(header)))  # 前7列是分类，第8列起是样品

    micro = [r for r in rows if r[0] in MICRO_KINGDOMS]
    macro = [r for r in rows if r[0] not in MICRO_KINGDOMS]
    print(f'总行数 {len(rows)} -> 微生物 {len(micro)} 行, 宏观生物 {len(macro)} 行')

    # 各自归一化（每个样品列除以该列小计）
    for table in (micro, macro):
        for si in samples:
            total = sum(float(r[si]) for r in table if r[si].strip())
            if total > 0:
                for r in table:
                    if r[si].strip():
                        r[si] = str(round(float(r[si]) / total, 10))

    micro_path = os.path.join(args.outdir, 'Abundance.micro.xls')
    macro_path = os.path.join(args.outdir, 'Abundance.macro.xls')
    write_table(micro_path, header, micro, samples)
    write_table(macro_path, header, macro, samples)

    if not args.keep_total:
        os.remove(args.input)
        print(f'已删除总表 {args.input}（如需保留请加 --keep-total）')
    else:
        print(f'总表已保留: {args.input}')
    print('完成: 输出 Abundance.micro.xls / Abundance.macro.xls')

if __name__ == '__main__':
    main()
