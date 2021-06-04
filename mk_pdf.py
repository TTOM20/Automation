import os
import img2pdf
from PIL import Image # img2pdfと一緒にインストールされたPillowを使います

if __name__ == '__main__':
    pdf_FileName = "./img/output.pdf" # 出力するPDFの名前
    png_Folder = "./img/" # 画像フォルダ
    extension  = ".png" # 拡張子がPNGのものを対象

    with open(pdf_FileName,"wb") as f:
        # 画像フォルダの中にあるPNGファイルを取得し配列に追加、バイナリ形式でファイルに書き込む
        f.write(img2pdf.convert([Image.open(png_Folder+j).filename for j in os.listdir(png_Folder)if j.endswith(extension)]))
