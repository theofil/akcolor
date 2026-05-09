#!/usr/bin/env python3
import pathlib
import ROOT

TDR_STYLE = pathlib.Path.home() / 'qplot' / 'setTDRStyle.C'
DATA_DIR   = pathlib.Path('.')
FIGS       = pathlib.Path.home() / 'www/files/akcolor'

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.LoadMacro(str(TDR_STYLE))
ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gStyle.SetPalette(ROOT.kBird)

FIGS.mkdir(parents=True, exist_ok=True)

HIST_CLASSES = {'TH1F', 'TH2F', 'TH2Poly'}

for root_file in sorted(DATA_DIR.glob('*slimfriend.root')):
    stem = root_file.stem
    print('Processing', root_file)

    f = ROOT.TFile(str(root_file))
    if not f or f.IsZombie():
        print('  Could not open', root_file)
        continue

    for key in f.GetListOfKeys():
        classname = key.GetClassName()
        if classname not in HIST_CLASSES:
            continue

        obj = key.ReadObj()
        hname = obj.GetName()
        obj.SetTitle('')
        is2d = classname in ('TH2F', 'TH2Poly')
        c = ROOT.TCanvas(f'c_{stem}_{hname}', '', 600 if is2d else 800, 600)
        if is2d:
            c.SetRightMargin(0.15)
            obj.Draw('COLZ')
            if classname == 'TH2Poly':
                obj.GetXaxis().SetTitle('#Delta y_{i}')
                obj.GetYaxis().SetTitle('#Delta#varphi_{i}')
                obj.GetXaxis().SetNdivisions(505)
                obj.GetYaxis().SetNdivisions(505)
        else:
            obj.Draw('HIST')
        outpath = str(FIGS / f'{hname}_{stem}.pdf')
        c.SaveAs(outpath)
        c.Close()
        print('  Saved', outpath)

    f.Close()

# update index.php
index_path = FIGS / 'index.php'
php = '''<?php
$files = glob("*.pdf");
sort($files);
echo "<html><body>";
echo "<h2>akcolor plots</h2><ul>";
foreach ($files as $f) {
    echo "<li><a href=\\"$f\\">$f</a></li>";
}
echo "</ul></body></html>";
?>
'''
index_path.write_text(php)
print('Updated', index_path)
