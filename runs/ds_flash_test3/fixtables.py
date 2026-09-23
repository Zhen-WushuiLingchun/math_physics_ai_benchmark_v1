import re, io

def is_cjk(ch):
    o = ord(ch)
    return (0x4E00 <= o <= 0x9FFF) or (0x3000 <= o <= 0x303F) or (0xFF00 <= o <= 0xFFEF) \
        or (0x2010 <= o <= 0x2027)

def wrap_cell(cell):
    lead = cell[:len(cell) - len(cell.lstrip())]
    body = cell.strip()
    if not body:
        return cell
    if not any(is_cjk(ch) for ch in body):
        return cell
    if body.startswith('{') and body.endswith('}'):
        return cell
    m = re.match(r'^(.*?)(\s*\\\\\s*)$', body, re.S)
    if m:
        body = '{' + m.group(1) + '}' + m.group(2)
    else:
        body = '{' + body + '}'
    return lead + body

for fn in ['sec2_A.tex', 'sec3_B.tex', 'sec4_C.tex', 'sec5_num.tex', 'solution.tex']:
    lines = io.open(fn, encoding='utf-8').read().split('\n')
    out = []
    intab = False
    for ln in lines:
        if '\\begin{tabular}' in ln:
            intab = True
        if intab and ('&' in ln):
            cells = re.split(r'(?<!\\)&', ln)
            cells = [wrap_cell(c) for c in cells]
            ln = '&'.join(cells)
        if '\\end{tabular}' in ln or '\\bottomrule' in ln:
            intab = False
        out.append(ln)
    io.open(fn, 'w', encoding='utf-8').write('\n'.join(out))
    print('done', fn)
