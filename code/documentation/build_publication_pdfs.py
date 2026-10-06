"""Compose publication PDFs from the supplied Markdown; no scientific run.

Developed with ChatGPT/Codex assistance (OpenAI), under Patrick Royer's direction.
See ../../licenses/DEL_v1_1_FR.md for rights on this original script.
Dependencies: reportlab, matplotlib, Pillow; Python standard library.
"""
from pathlib import Path
import argparse
import hashlib
import html
import io
import json
import re
import tempfile

import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager, mathtext
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                               Table, TableStyle, PageBreak, KeepTogether)

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['savefig.transparent'] = True
W, H = A4
MARGIN = 48
WIDTH = W - 2*MARGIN
INK = '#20364D'
BLUE = '#1F5B75'


def register_fonts():
    faces = [
        ('DYN-Regular','DejaVu Sans','normal','normal'),
        ('DYN-Bold','DejaVu Sans','bold','normal'),
        ('DYN-Italic','DejaVu Sans','normal','oblique'),
        ('DYN-BoldItalic','DejaVu Sans','bold','oblique'),
        ('DYN-Mono','DejaVu Sans Mono','normal','normal'),
    ]
    for name,family,weight,style in faces:
        prop=font_manager.FontProperties(family=family,weight=weight,style=style)
        pdfmetrics.registerFont(TTFont(name,font_manager.findfont(prop)))
    pdfmetrics.registerFontFamily('DYN-Regular',normal='DYN-Regular',bold='DYN-Bold',italic='DYN-Italic',boldItalic='DYN-BoldItalic')


def normalize_math(expr):
    expr=' '.join(expr.split())
    expr=re.sub(r'\\ge(?![A-Za-z])',r'\\geq',expr)
    expr=re.sub(r'\\le(?![A-Za-z])',r'\\leq',expr)
    expr=expr.replace(r'\pmod{\mathbb{Z}^2}',r'\ (\mathrm{mod}\ \mathbb{Z}^2)')
    expr=expr.replace(r'\widetilde P',r'\widetilde{P}')
    return expr


class MathCache:
    def __init__(self,path):
        self.path=path
        path.mkdir(parents=True,exist_ok=True)
        self.display_count=0
        self.inline_count=0

    def render(self,expr,size,color=INK):
        expr=normalize_math(expr)
        key=hashlib.sha256(('transparent-v1|'+color+'|'+str(size)+'|'+expr).encode()).hexdigest()
        target=self.path/(key+'.png')
        dpi=300
        if not target.exists():
            mathtext.math_to_image('$'+expr+'$',target,
                prop=font_manager.FontProperties(size=size),dpi=dpi,format='png',color=color)
        with PILImage.open(target) as im:
            width,height=im.size
        return target,width*72/dpi,height*72/dpi

    def inline(self,expr,size=10.1,color=INK):
        self.inline_count+=1
        p,w,h=self.render(expr,size,color)
        scale=min(1,14.0/h)
        return f'<img src="{p}" width="{w*scale:.3f}" height="{h*scale:.3f}" valign="middle"/>'

    def display(self,expr):
        self.display_count+=1
        p,w,h=self.render(expr,17)
        scale=min(1,(WIDTH-26)/w)
        img=Image(str(p),width=w*scale,height=h*scale)
        img.hAlign='LEFT'
        table=Table([[img]],colWidths=[WIDTH])
        table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,-1),colors.HexColor('#F3F7FA')),
            ('BOX',(0,0),(-1,-1),0.4,colors.HexColor('#C6D4DE')),
            ('LEFTPADDING',(0,0),(-1,-1),13),('RIGHTPADDING',(0,0),(-1,-1),13),
            ('TOPPADDING',(0,0),(-1,-1),10),('BOTTOMPADDING',(0,0),(-1,-1),10),
        ]))
        return KeepTogether([Spacer(1,3),table,Spacer(1,9)])


def styles():
    base=dict(fontName='DYN-Regular',fontSize=10.1,leading=15.4,textColor=colors.HexColor('#202A33'),spaceAfter=7.2,splitLongWords=False,allowWidows=False,allowOrphans=False)
    return {
        'body':ParagraphStyle('body',**base),
        'title':ParagraphStyle('title',**{**base,'fontName':'DYN-Bold','fontSize':20.5,'leading':26,'textColor':colors.HexColor(INK),'spaceBefore':7,'spaceAfter':12,'keepWithNext':True}),
        'h2':ParagraphStyle('h2',**{**base,'fontName':'DYN-Bold','fontSize':14.0,'leading':18,'textColor':colors.HexColor(INK),'spaceBefore':12,'spaceAfter':7,'keepWithNext':True}),
        'h3':ParagraphStyle('h3',**{**base,'fontName':'DYN-Bold','fontSize':11.5,'leading':16,'textColor':colors.HexColor(BLUE),'spaceBefore':9,'spaceAfter':5,'keepWithNext':True}),
        'cell':ParagraphStyle('cell',**{**base,'fontSize':8.7,'leading':13.5,'spaceAfter':0}),
        'cellhead':ParagraphStyle('cellhead',**{**base,'fontName':'DYN-Bold','fontSize':8.7,'leading':13.5,'spaceAfter':0,'textColor':colors.white}),
        'bullet':ParagraphStyle('bullet',**{**base,'leftIndent':13,'firstLineIndent':-10,'spaceAfter':5}),
        'code':ParagraphStyle('code',fontName='DYN-Mono',fontSize=8,leading=12,textColor=colors.HexColor('#263748'),spaceAfter=5,splitLongWords=False),
    }


def inline(text,cache,size=10.1,color=INK):
    tokens=[]
    def token(value):
        tokens.append(value)
        return f'DYNTOKEN{len(tokens)-1}END'
    text=re.sub(r'\$([^$]+)\$',lambda m:token(cache.inline(m.group(1),size,color)),text)
    text=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',lambda m:token(f'<link href="{html.escape(m.group(2),quote=True)}" color="{BLUE}">{html.escape(m.group(1))}</link>'),text)
    text=html.escape(text)
    text=re.sub(r' +([;:!?])',r'&#160;\1',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<b>\1</b>',text)
    text=re.sub(r'`([^`]+)`',r'<font name="DYN-Mono">\1</font>',text)
    for i,value in enumerate(tokens): text=text.replace(f'DYNTOKEN{i}END',value)
    return text


def table_widths(n,header):
    if n==5 and 'Modèle' in header: return [103,164,75,90,WIDTH-432]
    if n==4: return [72,126,92,WIDTH-290]
    if n==3: return [190,(WIDTH-190)/2,(WIDTH-190)/2]
    return [WIDTH/n]*n


def parse_markdown(text,cache):
    ss=styles(); story=[]; lines=text.splitlines(); i=0
    while i<len(lines):
        line=lines[i]
        if not line.strip(): i+=1; continue
        if line.strip()=='<!-- pagebreak -->': story.append(PageBreak()); i+=1; continue
        if line.strip()=='$$':
            i+=1; expr=[]
            while i<len(lines) and lines[i].strip()!='$$': expr.append(lines[i]); i+=1
            if i==len(lines): raise ValueError('Unclosed math block')
            story.append(cache.display('\n'.join(expr))); i+=1; continue
        if line.lstrip().startswith('#'):
            depth=len(line)-len(line.lstrip('#'))
            heading=line[depth:].strip()
            style='title' if depth==1 else ('h2' if depth==2 else 'h3')
            story.append(Paragraph(inline(heading,cache,ss[style].fontSize),ss[style])); i+=1; continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                values=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?',x) for x in values): rows.append(values)
                i+=1
            n=len(rows[0]); assert all(len(r)==n for r in rows)
            data=[[Paragraph(inline(x,cache,8.7,'white' if j==0 else INK),ss['cellhead' if j==0 else 'cell']) for x in r] for j,r in enumerate(rows)]
            t=Table(data,colWidths=table_widths(n,' '.join(rows[0])),repeatRows=1,hAlign='LEFT')
            t.setStyle(TableStyle([
                ('BACKGROUND',(0,0),(-1,0),colors.HexColor(INK)),
                ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F3F6F8')]),
                ('GRID',(0,0),(-1,-1),0.3,colors.HexColor('#AFBEC9')),
                ('VALIGN',(0,0),(-1,-1),'TOP'),
                ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
                ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
            ]))
            story.extend([t,Spacer(1,9)]); continue
        if line.startswith('    ') or line.startswith('```'):
            if line.startswith('```'):
                i+=1; code=[]
                while i<len(lines) and not lines[i].startswith('```'): code.append(lines[i]); i+=1
                i+=1
            else:
                code=[]
                while i<len(lines) and lines[i].startswith('    '): code.append(lines[i][4:]); i+=1
            for command in code:
                size=min(8.0,(WIDTH-20)/max(pdfmetrics.stringWidth(command,'DYN-Mono',1),1))
                style=ParagraphStyle('fitted_code',parent=ss['code'],fontSize=size,leading=12)
                p=Paragraph(html.escape(command).replace(' ','&#160;'),style)
                t=Table([[p]],colWidths=[WIDTH])
                t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),colors.HexColor('#F3F6F8')),('BOX',(0,0),(-1,-1),0.3,colors.HexColor('#C6D4DE')),('LEFTPADDING',(0,0),(-1,-1),10),('RIGHTPADDING',(0,0),(-1,-1),10),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
                story.extend([t,Spacer(1,5)])
            continue
        is_bullet=line.startswith('- ') or re.match(r'^\d+\. ',line)
        para=[line]; i+=1
        while i<len(lines) and lines[i].strip() and not (lines[i].startswith(('#','|','- ','    ','```')) or lines[i].strip() in ('$$','<!-- pagebreak -->') or re.match(r'^\d+\. ',lines[i])):
            para.append(lines[i]); i+=1
        joined=' '.join(x.strip() for x in para)
        if line.startswith('- '): joined='• '+joined[2:]
        story.append(Paragraph(inline(joined,cache),ss['bullet' if is_bullet else 'body']))
    return story


def build(path,text,title,footer,cache):
    doc=SimpleDocTemplate(str(path),pagesize=A4,leftMargin=MARGIN,rightMargin=MARGIN,
        topMargin=47,bottomMargin=47,title=title,author='Patrick Royer',
        subject='Dynagénèse - recherche indépendante assistée par IA',
        creator='ChatGPT/Codex-assisted documentary integration - ReportLab')
    pages=[]
    def frame(canvas,doc):
        pages.append(doc.page)
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor('#B8C7D1')); canvas.setLineWidth(.45)
        canvas.line(MARGIN,H-34,W-MARGIN,H-34)
        canvas.line(MARGIN,34,W-MARGIN,34)
        canvas.setFillColor(colors.HexColor('#647586')); canvas.setFont('DYN-Regular',7.2)
        canvas.drawString(MARGIN,H-27,'Patrick Royer | Dynagénèse | recherche indépendante')
        canvas.drawString(MARGIN,23,footer)
        canvas.drawRightString(W-MARGIN,23,str(doc.page))
        canvas.restoreState()
    doc.build(parse_markdown(text,cache),onFirstPage=frame,onLaterPages=frame)
    return len(pages)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--cache',type=Path)
    args=ap.parse_args(); root=args.root.resolve()
    cache_root=args.cache or Path(tempfile.mkdtemp(prefix='dyn_pdf_'))
    register_fonts(); cache=MathCache(cache_root)
    main_md=root/'manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.md'
    main_pdf=main_md.with_suffix('.pdf')
    pages=build(main_pdf,main_md.read_text(),
        'Déficit de périmètre des enveloppes entières randomisées',
        'v0.1b | conjecture à deux termes | assistance IA | DEL 1.1',cache)
    display_count=cache.display_count
    fr=(root/'licenses/DEL_v1_1_FR.md').read_text()
    en=(root/'licenses/DEL_v1_1_EN.md').read_text()
    licence_pdf=root/'licenses/DYN-LICENCE_DEL_v1_1.pdf'
    licence_pages=build(licence_pdf,fr+'\n\n<!-- pagebreak -->\n\n'+en,
        'Licence Dynagénèse Éthique (DEL), version 1.1',
        'DEL 1.1 | texte français de référence | traduction anglaise de lecture',cache)
    result={'manuscript_pages':pages,'display_equations':display_count,'licence_pages':licence_pages,
            'inline_math_occurrences':cache.inline_count,
            'manuscript_sha256':hashlib.sha256(main_pdf.read_bytes()).hexdigest(),
            'licence_sha256':hashlib.sha256(licence_pdf.read_bytes()).hexdigest(),
            'scope':'Documentary rendering only; no scientific generator or P2.'}
    (root/'evidence/PDF_BUILD_v0_1b.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
