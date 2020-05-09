#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import re
import os
import sys
import platform
import logging
import hashlib
import html_util
from io import open

# initialize the logger
logger = logging.getLogger('generate_rem_report')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

filterQCList= ["Reads_with_N", "Reads_with_N_Rate", "Reads_with_LowQuality", "Reads_with_LowQuality_Rate",
    "Reads_with_Adapter", "Reads_with_Adapter_Rate", "Reads_with_Duplications", "Reads_with_Duplications_Rate",
    "Raw_Base_Q20", "Clean_Base_Q20"]

def getFieldNameByLang(fieldName, lang="cn"):
    if lang == "cn":
        fieldDict = {
            "Multi-PCR SARS-CoV-2 Analysis Report": "多重PCR新冠病毒分析报告",
            "English": "中文",
            "Basic summary": "基本信息",
            "QC result": "质控结果",
            "Result": "表格链接",
            "Sample": "样本",
            "Identification and Quantification": "鉴定结果",
            "Total_Reads": "Reads总数",
            "No.": "序号",
            "Reads_Number": "reads数",
			"SARS-CoV-2":"新型冠状病毒",
			"Assembly":"组装",
            "Depth Distribution":"深度分布",
            "GC_Content":"GC含量",
            "Raw_Reads":"原始Reads",
            "Mapping_Rate":"比对率",
            "SARS-CoV-2_Reads_Number":"SARS-CoV-2 Reads数目",
            "SARS-CoV-2_Reads_Pct":"SARS-CoV-2比例",
            "≥1X_Cov":"≥1X覆盖度",
            "≥100X_Cov":"≥100X覆盖度",
            "Identification_Result":"阴阳性",
            "Picture":"图片链接",
            "Consensus Sequence":"一致性序列",
            "VCF":"VCF文件",
			"Download":"下载",
			"Positive":"阳性",
			"Negative":"阴性",
            }
        if fieldName in fieldDict:
            fieldName = fieldDict[fieldName]
    return fieldName

def setWarpTitle(title):
    titleName = ""
    if "/" in title:
        titleName = title.split("/")[0] + "<br>"
        for name in title.split("/")[1:]:
            titleName += "/" + name + "<br>"
    elif "_" in title:
        for name in title.split("_"):
            titleName += name + "<br>"
    else:
        titleName = title
    return titleName

def getTableContent(fileName, isWarpTitle = True, lang = 'cn'):
    try:
        content = ""
        with open(fileName, 'r') as fh:
            rowNum = 0
            for line in fh:
                line = line.strip()
                rowNum += 1
                # filter the items that will not be displayed in the report
                if rowNum != 1 and line.split("\t")[0] in filterQCList:
                    continue
                
                content += "<tr>"
                if rowNum == 1:
                    for field in line.split("\t"):
                        field = getFieldNameByLang(field, lang)
                        if isWarpTitle:
                            field = setWarpTitle(field)
                        content += '<th scope="col">%s</th>' % field
                else:
                    for field in line.split("\t"):
                        if field.isdigit():
                            field = format(int(field), ',')
                        content += '<td scope="row">%s</td>' % field
                content += "</tr>"
    except:
        content = '<div class="noDataTitle">无数据文件: %s</div>' % fileName  
    return content

def generate_html_report(path, name='rem', output_path='E:/codezlims/rem/Result/sample1', lang="cn"):
    html=   '''
            <!DOCTYPE html>
            <html>
            <head>
                <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
                <title>'''+ getFieldNameByLang("Multi-PCR SARS-CoV-2 Analysis Report", lang) + '''</title>
                <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">
                <meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1, maximum-scale=1, user-scalable=no">
                <link href="https://fonts.googleapis.com/css?family=Roboto" rel="stylesheet">
                <!-- load css file -->
                <style type="text/css">
            ''' + html_util.HtmlUtil.getFileContent('base.css') + '''
            ''' + html_util.HtmlUtil.getFileContent('common.css') + '''
            ''' + html_util.HtmlUtil.getFileContent('table.css') + '''
                </style>
            </head>
            <body>
            <!--header starts-->
            <div class="header">
                <div class="repeater"></div>
                <div class="wrapper">
                    <div class="headerBox clearfloat">
                        <div class="headLeft fl">
                        <h1>'''+ getFieldNameByLang("Multi-PCR SARS-CoV-2 Analysis Report", lang) + '''</h1>
                        <div style="float:right; margin:-40px -130px 10px 160px; background-color: #f7faeb; width: 100px; height: 30px; border-radius: 8px;text-align:center;"><a href="./'''+ name 
                        
    if lang == "cn":
        html+= "_en"
    else:
        html+= "_cn"
    html+=               '''.html" style="color: #1c567f; font-size: 20px; padding-top: 2px;">'''
    if lang == "cn":
        html+= "English"
    else:
        html+= "中文"
    html+=          '''</a></div>
                        <h2>V1.0.0.0</h2>
                        </div>
                        <!-- headLeft -->
                        <div class="headRight fr">
                        <div class="logo">
                        <img src="''' + html_util.HtmlUtil.getPNGBinary('logo.png') + '''">
                        </div>
                        </div>
                    <!-- headRight -->
                    </div>
                </div>
            </div>
            <!--header ends-->

            <!--container starts-->
            <div class="container">
            <div class="wrapper">
            <div class="secOne">
                <h1 class="headText headTextOne">
                    1. '''+ getFieldNameByLang("Basic summary", lang) + '''
                    <div class="headTextIconDIV">
                        <img id="iDivOne" src="''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''" onclick="showAndHidden1();">
                    </div>
                </h1>
                <div id="divOne">
					<div class="secBox">
						<h1>''' + getFieldNameByLang("Sample", lang) + '''</h1>
					</div>
					<div class="content">
					<p>''' + name +'''</p></div>
                    <div class="secBox" style="border: none">
                        <h1>'''+ getFieldNameByLang("QC result", lang) + '''</h1>
                    </div>
                    <div class="content">
                    <p>'''
    if lang == "cn":
        html+=  '''对原始测序reads进行质控，去除含有接头的、低质量碱基数目比例超过阈值的、带有‘N’碱基数目比例超过阈值的reads，最后统计出Clean Reads及其相应的比例。'''
    else: 
        html+=  '''We remove the low-quality reads from raw FASTQ by removing adapter contained sequences, N base ratio and low-quality base exceed the thresholds.'''    
    html+=      '''</p>
                    </div>
                    <div class="dataLink">
                        <a href="./QC.xlsx" download>'''+ getFieldNameByLang("Result", lang) + '''</a>
                    </div>
                    <div class="secThree">
                        <table cellspacing="0" style="width: 510px;">
                            <tbody>'''
    html += getTableContent(os.path.join(path, 'QC.txt'), False, lang)  
    html +=             '''</tbody>
                        </table>
                    </div>
                    <div class="content"><p></p></div>
                    <div class="content">'''
    if lang =="cn":
        html += '''     
                        <p>1. 样本：样本名；</p>
                        <p>2. Raw_Q30: 原始Reads的Q30；</p>
                        <p>3. GC_Content：原始数据GC含量；</p>
                        <p>4. 原始Reads：原始reads数；</p>
                        <p>5. Clean_Reads：clean reads数；</p>
                        <p>6. Clean_Rate：clean reads数/原始reads数*100%；</p>
                        <p>7. 比对率：与3个目标参考序列（lambda DNA，GAPDH和新冠病毒）总的比对率；</p>
                '''
    else:
        html += '''
                        <p>1. Sample: Sample name;</p>
                        <p>2. Raw_Q30: The ratio of base whose quality > 30 of the raw FASTQ;</p>
                        <p>3. GC_Content：The GC bases proportion of raw FASTQ;</p>
                        <p>4. Raw_Reads：The number of raw FASTQ;</p>
                        <p>5. Clean_Reads：The number of clean FASTQ;</p>
                        <p>6. Clean_Rate：The ratio of clean FASTQ/raw FASTQ;</p>
                        <p>7. Mapping_Rate: The mapping ratio of clean fastq maps against three reference sequences (lambda DNA, GAPDH, coronavirus);</p>
                '''

    html += ''' 
                </div>
                </div>             
                <script type="text/javascript">
                var div = document.getElementById('divOne');
                div.style.display = 'block';

                function showAndHidden1() {
                    if (div.style.display == 'block') {
                    div.style.display = 'none';
                    document.getElementById('iDivOne').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-down.png') + '''";

                    } else {
                    div.style.display = 'block';
                    document.getElementById('iDivOne').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''";
                    }
                }
                </script>
            </div>
            <!-- secTwoMain -->
            <div class="secOne">
                <h1 class="headText headTextOne">2. '''+ getFieldNameByLang("Identification and Quantification", lang) + '''
                    <div class="headTextIconDIV"><img id="iDivTwo" src="''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''" onclick="showAndHidden2();"></div>
                </h1>
                <div id="divTwo">
                    <div class="secBox">
                        <h1>'''+ getFieldNameByLang("SARS-CoV-2", lang) + '''</h1>
                    </div>
                    <div class="content">
                        <p>'''
    if lang == "cn":
        html += '''使用BWA进行比对，使用SAMtools进行比对统计。'''
    else:
        html += '''BWA is used for mapping and SAMtools is used for mapping statistics.'''
    html +=            '''</p>
                    </div>
                    <div class="dataLink">
                        <a href="./Identification.xlsx" download>'''+ getFieldNameByLang("Result", lang) + '''</a>
                    </div>           
                <div class="secThree" style="margin-bottom:20px; ">
                    <table cellspacing="0" style="width: 1080px;">
                        <tbody>'''
    html += getTableContent(os.path.join(path, 'Identification.txt'), True, lang)                     
    html +=        '''</tbody>
                    </table>
                </div>
                <div class="content">
                    '''
    if lang == "cn":                
        html += '''<p>1. 样本：样品名；</p>
                   <p>2. Clean_Reads：clean reads数；</p>
                   <p>3. SARS-CoV-2_Reads：鉴定出的新冠病毒read数；</p>
                   <p>4. SARS-CoV-2_Reads_Pct：鉴定出的新冠病毒比例，等于100%*(SARS-CoV-2 reads number)/((lambda DNA reads number)+(SARS-CoV-2 reads number))；</p>
                   <p>5. ≥1X_Coverage：鉴定出的新冠病毒reads与reference的比对后统计≥1X的覆盖度；</p>
                   <p>6. ≥100X_Coverage：鉴定出的新冠病毒reads与reference的比对后统计≥100X的覆盖度；</p>
                   <p>7. Identification_Result：SARS-CoV-2_Reads_Pct大于或等于0.1%，为阳性；0.05%-0.1%为灰区；小于0.05为阴性；</p>
                '''
    else:
        html += '''<p>1. Sample: Sample name;</p>
                   <p>2. Clean_Reads: The number of clean FASTQ;</p>
                   <p>3. SARS-CoV-2_Reads: The reads number of identified SARS-CoV-2;</p>
                   <p>4. SARS-CoV-2_Reads_Pct: The proportion of SARS-CoV-2 reads: 100%*(SARS-CoV-2 reads number)/((lambda DNA reads number)+(SARS-CoV-2 reads number));</p>
                   <p>5. ≥1X_Coverage: The coverage is calculated as proportion of reference genome covered by more than 1 folds reads;</p>
                   <p>6. ≥100X_Coverage: The coverage is calculated as proportion of reference genome covered by more than 100 folds reads;</p>
                   <p>7. Identification_Result: Sample identification result according to SARS-CoV-2_Reads_Pct: Positive(≥0.1%), Indetermination(0.05%-0.1%), Negative(<0.05%);</p>
                '''  
    html +=     '''
                </div>   
                </div>
            <script type="text/javascript">
                var divTwo = document.getElementById('divTwo');
                divTwo.style.display = 'block';

                function showAndHidden2() {
                    if (divTwo.style.display == 'block') {
                    divTwo.style.display = 'none';
                    document.getElementById('iDivTwo').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-down.png') + '''";

                    } else {
                    divTwo.style.display = 'block';
                    document.getElementById('iDivTwo').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''";
                    }
                }
            </script>
        </div>'''
    svg_file = open(os.path.join(path, 'Windows.Depth.svg'), 'r')
    svg_string = svg_file.read()
    html += '''
        <!-- secThreeMain -->
        <div class="secOne">
            <h1 class="headText headTextOne">
                3. ''' + getFieldNameByLang("Assembly", lang) + '''
                <div class="headTextIconDIV">
                    <img id="iDivThree" src="''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''" onclick="showAndHidden3();">
                </div>
            </h1>
            <div id="divThree">
               <div class="secBox">
                   <h1>''' + getFieldNameByLang("Consensus Sequence", lang) + '''</h1>
               </div>
               <div class="dataLink">
                   <a href="./''' + name + '''.Consensus.fa" download>'''+ getFieldNameByLang("Download", lang) + '''</a>
               </div>
               <div class="content">
                   <p>
               '''
    if lang == "cn":
        html += '''使用bcftools构建一致性序列。'''
    else:
        html += '''BCFtools is used to create consensus sequence by applying VCF variants.'''
    html += '''
                   </p>
               </div>
               <div class="secBox" style="border: none">
                   <h1>''' + getFieldNameByLang("VCF", lang) + '''</h1>
               </div>
               <div class="dataLink">
                   <a href="./''' + name + '''.vcf.gz" download>'''+ getFieldNameByLang("Download", lang) + '''</a>
               </div>
               <div class="content">
                   <p>
               '''
    if lang == "cn":
        html += '''使用freebayes进行变异检测。'''
    else:
        html += '''freebayes is used for find small polymorphisms.'''
    html += '''
                   </p>
               </div>
               <div class="secBox" style="border: none">
                   <h1>''' + getFieldNameByLang("Depth Distribution", lang) + '''</h1>
               </div>
               <div class="dataLink">
                   <a href="./Windows.Depth.svg" download>'''+ getFieldNameByLang("Download", lang) + '''</a>
               </div>
               <div class="svgPic">''' + svg_string + '''
               </div>
            </div>
        <script type="text/javascript">
            var divThree = document.getElementById('divThree');
            divThree.style.display = 'block';

            function showAndHidden3() {
                if (divThree.style.display == 'block') {
                    divThree.style.display = 'none';
                    document.getElementById('iDivThree').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-down.png') + '''";
                } else {
                    divThree.style.display = 'block';
                    document.getElementById('iDivThree').src = "''' + html_util.HtmlUtil.getPNGBinary('arrow-up.png') + '''";
                }
            }
        </script>
        </div>
    '''
    html += '''
        </div>
        <div class="repeater"></div>
        <div class="secBottom"></div>
        </div>
        </body>
        </html>
    '''
    with open(output_path + '/' + name + '_' + lang +'.html', 'w+', encoding='utf-8') as report:
        print(output_path + '/' + name + '_' + lang +'.html') 
        report.write(html)


if __name__ == '__main__':

    if len(sys.argv) < 2:
        logger.error("Please provide the path of the result folder.")
        sys.exit(-1)
    path = sys.argv[1]
#    name = os.path.basename(path)

#    if not name:
#        name = os.path.basename(path[:-1])
    name = sys.argv[2]

    logger.info('\nResult folder is: %s\nSample name is: %s' % (path, name))
    # Local test
    generate_html_report(path, name, path, "en")
    generate_html_report(path, name, path, "cn")
