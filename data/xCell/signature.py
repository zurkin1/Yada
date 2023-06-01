ctypes = ['Monocytes', 'CD8+ T-cells', 'NK cells', 'Macrophages', 'Endothelial cells', 'Neutrophils', 'Erythrocytes', 'CD8+ naive T-cells', 'CD4+ naive T-cells', 'MPP', 'Smooth muscle',	
          'Fibroblasts', 'Epithelial cells', 'Epithelial cells', 'Keratinocytes', 'Chondrocytes', 'Adipocytes', 'B-cells', 'CD4+ T-cells', 'CD8+ Tem', 'CMP', 'GMP', 'Tregs', 'HSC',
          'Plasma cells', 'B-cells', 'CD4+ Tcm', 'mv Endothelial cells', 'CD4+ Tem', 'Memory B-cells', 'CD8+ Tcm', 'naive B-cells', 'Eosinophils', 'Macrophages M1', 'Myocytes', 
          'ly Endothelial cells', 'MSC', 'Macrophages M2', 'Osteoblasts', 'Preadipocytes', 'Melanocytes', 'Skeletal muscle', 'CD4+ memory T-cells', 'Megakaryocytes', 
          'pro B-cells', 'Basophils', 'cDC', 'Astrocytes', 'Neurons', 'NKT', 'Osteoblast', 'Pericytes', 'Platelets', 'Mast cells', 'Hepatocytes', 'MEP', 'Mesangial cells', 'CLP',
          'pDC', 'aDC', 'DC', 'Sebocyres', 'Tgd cells', 'Th1 cells', 'Th2 cells', 'Sebocytes']

newcols = [0]*len(cols)
i = 0
for col in cols:
   for cell in ctypes:
       if cell in col:
          newcols[i] = cell
          continue
   i = i + 1