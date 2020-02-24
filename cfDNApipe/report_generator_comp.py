# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 15:54:49 2019
@author: Jiaqi Huang
"""

from yattag import Doc, indent
from .Configure import Configure
from .Configure2 import Configure2
import datetime
import os
import shutil
import bz2
import pkg_resources


def report_generator_comp(
    case_fastqcRes=None,
    case_identifyAdapterRes=None,
    case_bismarkRes=None,
    case_deduplicateRes=None,
    case_rmduplicateRes=None,
    case_fraglenplotRes=None,
    ctrl_fastqcRes=None,
    ctrl_identifyAdapterRes=None,
    ctrl_bismarkRes=None,
    ctrl_deduplicateRes=None,
    ctrl_rmduplicateRes=None,
    ctrl_fraglenplotRes=None,
    OCFRes=None,
    CNVRes=None,
    outputdir=None,
    label=None, #list
):

    if outputdir is None:
        outputdir = Configure2.getRepDir()
    if label is None:
        label = ["Case", "Control"]

    doc, tag, text, line = Doc().ttl()
    write_head(doc, tag, text, line)
    write_body(
        doc,
        tag,
        text,
        line,
        case_fastqcRes,
        case_identifyAdapterRes,
        case_bismarkRes,
        case_deduplicateRes,
        case_rmduplicateRes,
        case_fraglenplotRes,
        ctrl_fastqcRes,
        ctrl_identifyAdapterRes,
        ctrl_bismarkRes,
        ctrl_deduplicateRes,
        ctrl_rmduplicateRes,
        ctrl_fraglenplotRes,
        OCFRes,
        CNVRes,
        outputdir,
        label,
    )

    fout = open(os.path.join(outputdir, "Cell Free DNA WGBS Analysis Report.html"), "w")
    fout.write(indent(doc.getvalue()))
    fout.close()


def write_head(doc, tag, text, line):
    # read the header bz2 file
    href_file = pkg_resources.resource_filename("cfDNApipe", "data/src_href.bz2")
    with bz2.open(href_file, "rt") as fscript:
        textscript = fscript.read()
    srch = textscript.split("\n")

    # write the html head
    doc.asis("<!DOCTYPE html>")

    with tag("html", xmlns="http://www.w3.org/1999/xhtml"):
        with tag("head"):
            doc.stag("meta", charset="utf-8")
            doc.asis(
                '<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />'
            )
            doc.stag("meta", name="date", content=str(datetime.date.today()))

            with tag("title"):
                text("Cell Free DNA WGBS Analysis Report")

            with tag("script", src=srch[0]):
                text()

            doc.stag(
                "meta", name="viewport", content="width=device-width, initial-scale=1"
            )

            with tag("link", href=srch[1], rel="stylesheet"):
                text()

            with tag("script", src=srch[2]):
                text()

            with tag("script", src=srch[3]):
                text()

            with tag("script", src=srch[4]):
                text()

            with tag("script", src=srch[5]):
                text()

            with tag("link", href=srch[6], rel="stylesheet"):
                text()

            with tag("script", src=srch[7]):
                text()

            with tag("script", src=srch[8]):
                text()

            with tag("link", href=srch[9], rel="stylesheet"):
                text()

            with tag("script", src=srch[10]):
                text()

            with tag("link", href=srch[11], rel="stylesheet"):
                text()

            with tag("script", src=srch[12]):
                text()

            with tag("style", type="text/css"):
                text("code{white-space: pre;}")

            with tag("style", type="text/css"):
                text("\n  pre:not([class]){\n    background-color: white;\n  }")

            with tag("script", type="text/javascript"):
                doc.asis(
                    '\nif (window.hljs) {\n  hljs.configure({languages: []});\n  hljs.initHighlightingOnLoad();\n  if (document.readyState && document.readyState === "complete") {\n    window.setTimeout(function() { hljs.initHighlighting(); }, 0);\n  }\n}'
                )

            with tag("style", type="text/css"):
                text("\nh1 {\n  font size: 34px;\n}\n")
                text("h1.title {\n  font size: 38px;\n}\n")
                text("h2 {\n  font size: 30px;\n}\n")
                text("h3 {\n  font size: 24px;\n}\n")
                text("h4 {\n  font size: 18px;\n}\n")
                text("h5 {\n  font size: 16px;\n}\n")
                text("h6 {\n  font size: 12px;\n}\n")
                text(".table th:not([align]) {\n  text-align: left;\n}")

            with tag("style", type="text/css"):
                text(
                    "\ntable.customize {\n	font-size:12px;\n	color:#333333;\n	border-width: 1px;\n	border-color: #888888;\n	border-collapse: collapse;\n}"
                )
                text(
                    "\ntable.customize th {\n	border-width: 1px;\n	padding: 8px;\n	border-style: solid;\n	border-color: #888888;\n}"
                )
                text(
                    "\ntable.customize td {\n	border-width: 2px;\n	padding: 6px;\n	border-left-style: none;\n	border-right-style: none;\n	border-top-style: solid;\n	border-bottom-style: solid;\n	border-color: #DFDFDF;\n}"
                )


def write_body(
    doc,
    tag,
    text,
    line,
    case_fastqcRes,
    case_identifyAdapterRes,
    case_bismarkRes,
    case_deduplicateRes,
    case_rmduplicateRes,
    case_fraglenplotRes,
    ctrl_fastqcRes,
    ctrl_identifyAdapterRes,
    ctrl_bismarkRes,
    ctrl_deduplicateRes,
    ctrl_rmduplicateRes,
    ctrl_fraglenplotRes,
    OCFRes,
    CNVRes,
    outputdir,
    label,
):
    with tag("body"):
        with tag("style", type="text/css"):
            text(
                "\n.main-container {\n  max-width: 940px;\n  margin-left: auto;\n  margin-right: auto;\n}"
            )
            text(
                "\ncode {\n  color: inferit;\n  background-color: rgba(0, 0, 0, 0.04);\n}"
            )
            text("\nimg {\n  max-width: 100%;\n  height: auto;\n}")
            text("\n.tabbed-pane {\n  padding-top: 12px;\n}")
            text("\nbutton.code-folding-btn:focus {\n  outling: none;\n}")

        with tag("div", klass="container-fluid main-container"):
            with tag("script"):
                doc.asis(
                    '\n$(document).ready(function () {\n  window.buildTabsets("TOC");\n});'
                )

            with tag("script"):
                doc.asis(
                    "\n$(document).ready(function ()  {\n\n    // move toc-ignore selectors from section div to header\n    $('div.section.toc-ignore')\n        .removeClass('toc-ignore')\n        .children('h1,h2,h3,h4,h5').addClass('toc-ignore');\n\n    // establish options\n    var options = {\n      selectors: \"h1,h2,h3\",\n      theme: \"bootstrap3\",\n      context: '.toc-content',\n      hashGenerator: function (text) {\n        return text.replace(/[.\\/?&!#<>]/g, '').replace(/\s/g, '_').toLowerCase();\n      },\n      ignoreSelector: \".toc-ignore\",\n      scrollTo: 0\n    };\n    options.showAndHide = true;\n    options.smoothScroll = true;\n\n    // tocify\n    var toc = $(\"#TOC\").tocify(options).data(\"toc-tocify\");\n});"
                )

            with tag("style", type="text/css"):
                text(
                    "\n#TOC {\n  margin: 25px 0px 20px 0px;\n}\n@media (max-width: 768px) {\n#TOC {\n  position: relative;\n  width: 100%;\n}\n}"
                )
                text(
                    "\n.toc-content {\n  padding-left: 30px;\n  padding-right: 40px;\n}"
                )
                text("\ndiv.main-container {\n  max-width: 1200px;\n}")
                text(
                    "\ndiv.tocify {\n  width: 20%;\n  max-width: 260px;\n  max-height: 85%;\n}"
                )
                text(
                    "\n@media (min-width: 768px) and (max-width: 991px) {\n  div.tocify {\n    width: 25%;\n  }\n}"
                )
                text(
                    "\n@media (max-width: 767px) {\n  div.tocify {\n    width: 100%;\n    max-width: none;\n  }\n}"
                )
                text("\n.tocify ul, .tocify li {\n  line-height: 20px;\n}")
                text(
                    "\n.tocify-subheader .tocify-item {\n  font-size: 0.90em;\n  padding-left: 25px;\n  text-indent: 0;\n}"
                )
                text("\n.tocify .list-group-item {\n  border-radius: 0px;\n}")

            with tag("div", klass="row-fluid"):
                with tag("div", klass="col-xs-12 col-sm-4 col-md-3"):
                    with tag("div", id="TOC", klass="tocify"):
                        text()

                with tag("div", klass="toc-content col-xs-12 col-sm-8 col-md-9"):
                    with tag("div", klass="fluid-row", id="header"):
                        with tag("h1", klass="title toc-ignore"):
                            text("Cell Free DNA WGBS Analysis Report")

                        with tag("h4", klass="date"):
                            with tag("em"):
                                text(str(datetime.date.today()))

                #Case section
                    with tag(
                            "div",
                            id="Case",
                            klass="section level1",
                            style="margin:20px",
                        ):
                        with tag("h1"):
                            with tag("span", klass="header-section-number"):
                                text("1 " + label[0] + " analysis")
                        text("Subsections below are the results of the steps applied on the " + label[0] + " data.")
                        
                        case_title_count = 1
                        
                        # fastqc report
                        if case_fastqcRes is not None:
                            with tag(
                                "div",
                                id="case_fastqc_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Fastq Quality Control")
                                write_fastqc_report(
                                    doc, tag, text, line, case_fastqcRes, outputdir
                                )
                            case_title_count += 1

                        # identifyadapter report
                        if case_identifyAdapterRes is not None:
                            with tag(
                                "div",
                                id="case_idadapter_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Adapter Detection")
                                write_identifyadapter_report(
                                    doc, tag, text, line, case_identifyAdapterRes
                                )
                            case_title_count += 1

                        # bismark report
                        if case_bismarkRes is not None:
                            with tag(
                                "div",
                                id="case_bismark_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Bismark Alignment")
                                write_bismark_report(doc, tag, text, line, case_bismarkRes)
                            case_title_count += 1
                        
                        # deduplicate report
                        if case_deduplicateRes is not None:
                            with tag(
                                "div",
                                id="case_deduplicate_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Deduplicate Alignment")
                                write_deduplicate_report(doc, tag, text, line, case_deduplicateRes)
                            case_title_count += 1

                        # rmduplicate report
                        if case_rmduplicateRes is not None:
                            with tag(
                                "div",
                                id="case_rmduplicate_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Remove Duplicates")
                                write_rmduplicate_report(
                                    doc, tag, text, line, case_rmduplicateRes
                                )
                            case_title_count += 1

                        # fraglenplot report
                        if case_fraglenplotRes is not None:
                            with tag(
                                "div",
                                id="case_fraglenplot_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("1." + str(case_title_count))

                                    text(" " + label[0] + " Fragment Length Distribution")
                                write_fraglenplot_report(
                                    doc, tag, text, line, case_fraglenplotRes, outputdir
                                )
                            case_title_count += 1
                            
                #Control section
                    with tag(
                        "div",
                        id="Ctrl",
                        klass="section level1",
                        style="margin:20px",
                        ):
                        
                        with tag("h1"):
                            with tag("span", klass="header-section-number"):
                                text("2 " + label[1] + " analysis")
                        
                        text("Subsections below are the results of the steps applied on the " + label[1] + " data.")
                        
                        ctrl_title_count = 1
                        
                        # fastqc report
                        if ctrl_fastqcRes is not None:
                            with tag(
                                "div",
                                id="ctrl_fastqc_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("2." + str(ctrl_title_count))

                                    text(" " + label[1] + " Fastq Quality Control")
                                write_fastqc_report(
                                    doc, tag, text, line, ctrl_fastqcRes, outputdir
                                )
                            ctrl_title_count += 1

                        # identifyadapter report
                        if ctrl_identifyAdapterRes is not None:
                            with tag(
                                "div",
                                id="ctrl_idadapter_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("2." + str(ctrl_title_count))

                                    text(" " + label[1] + " Adapter Detection")
                                write_identifyadapter_report(
                                    doc, tag, text, line, ctrl_identifyAdapterRes
                                )
                            ctrl_title_count += 1

                        # bismark report
                        if ctrl_bismarkRes is not None:
                            with tag(
                                "div",
                                id="ctrl_bismark_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("2." + str(ctrl_title_count))

                                    text(" " + label[1] + " Bismark Alignment")
                                write_bismark_report(doc, tag, text, line, ctrl_bismarkRes)
                            ctrl_title_count += 1

                        # rmduplicate report
                        if ctrl_rmduplicateRes is not None:
                            with tag(
                                "div",
                                id="ctrl_rmduplicate_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("2." + str(ctrl_title_count))

                                    text(" " + label[1] + " Remove Duplicates")
                                write_rmduplicate_report(
                                    doc, tag, text, line, ctrl_rmduplicateRes
                                )
                            ctrl_title_count += 1

                        # fraglenplot report
                        if ctrl_fraglenplotRes is not None:
                            with tag(
                                "div",
                                id="ctrl_fraglenplot_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                with tag("h2"):
                                    with tag("span", klass="header-section-number"):
                                        text("2." + str(ctrl_title_count))

                                    text(" " + label[1] + " Fragment Length Distribution")
                                write_fraglenplot_report(
                                    doc, tag, text, line, ctrl_fraglenplotRes, outputdir
                                )
                            ctrl_title_count += 1

                #Compare section
                    comp_title_count = 3
                    
                    #OCF report
                    if OCFRes is not None:
                        with tag(
                            "div",
                            id="Comp-OCF",
                            klass="section level1",
                            style="margin:20px",
                            ):
                            with tag("h1"):
                                with tag("span", klass="header-section-number"):
                                    text(str(comp_title_count))
                                text(" OCF Results")
                            with tag(
                                "div",
                                id="OCF_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                write_OCF_report(
                                    doc, tag, text, line, OCFRes, outputdir
                                )
                        comp_title_count += 1
                    
                    #CNV report
                    if CNVRes is not None:
                        with tag(
                            "div",
                            id="Comp-CNV",
                            klass="section level1",
                            style="margin:20px",
                            ):
                            with tag("h1"):
                                with tag("span", klass="header-section-number"):
                                    text(str(comp_title_count))
                                text(" CNV Results")
                            with tag(
                                "div",
                                id="CNV_report",
                                klass="section level2",
                                style="margin:20px",
                            ):
                                write_CNV_report(
                                    doc, tag, text, line, CNVRes, outputdir, label
                                )
                        comp_title_count += 1
                        
        with tag("script"):
            doc.asis(
                "\nfunction bootstrapStylePandocTables() {\n  $('tr.header').parent('thead').parent('table').addClass('table table-condensed');\n$(document).ready(function () {\n  bootstrapStylePandocTables();\n});"
            )

        with tag("script"):
            doc.asis(
                '\n (function () {\n    var script = document.createElement("script");\n    script.type = "text/javascript";\n    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";\n    document.getElementsByTagName("head")[0].appendChild(script);\n  })();'
            )

def write_bismark_report(doc, tag, text, line, report_dir, max_sample=3):
    sample_num = 0
    for report in report_dir.getOutput("bismkRepOutput"):
        sample_num += 1
        if sample_num > max_sample:  # ignore the rest to shorten the report length
            break
        with tag("div", id="bismark_report_sub", klass="section level3"):
            with tag("h3"):
                text("Sample: " + report.split("/")[-1].split(".")[0])
            write_bismark_report_contents(doc, tag, text, line, report)

def write_bismark_report_contents(doc, tag, text, line, report):
    # read and locate PE_report.txt
    fin = open(report, "r", encoding="utf-8")
    cont = fin.readlines()
    for i in range(len(cont)):
        cont[i] = cont[i].strip("\n")
        if cont[i] == "Final Alignment report":
            p = i

    # contents
    with tag("div", style="line-height:30px"):
        with tag("ul", id="bismark_report_content"):
            for j in range(p + 2, p + 8):
                line("li", cont[j])

            with tag("li"):
                text(cont[p + 9])

                with tag("div", style="line-height:20px"):
                    with tag("table", klass="customize", width="100%"):
                        for j in range(p + 10, p + 14):
                            linex = cont[j].split("\t")
                            with tag("tr"):
                                if j % 2 == 0:
                                    doc.attr(klass="odd")
                                else:
                                    doc.attr(klass="even")

                                with tag("td", align="left"):
                                    text(linex[0])

                                with tag("td", align="left"):
                                    text(linex[1])

                                with tag("td", align="left"):
                                    text(linex[2])

                line("li", cont[p + 15])

    fin.close()

def write_identifyadapter_report(doc, tag, text, line, report_dir, max_sample=3):
    sample_num = 0
    for tmp_output in report_dir.getOutputs():
        if "-adapterFile" in tmp_output:
            sample_num += 1
            if sample_num > max_sample:  # ignore the rest to shorten the report length
                break
            report = report_dir.getOutput(tmp_output)
            with tag("div", id="idadapters_sub", klass="section level3"):
                with tag("h3"):
                    text(
                        "Sample: " + report.split("/")[-1].replace("-adapters.log", "")
                    )
                write_identifyadapter_report_contents(doc, tag, text, line, report)

def write_identifyadapter_report_contents(doc, tag, text, line, report):
    # read and locate adapters.log
    fin = open(report, "r", encoding="utf-8")
    cont = fin.readlines()

    # contents
    with tag("div", style="line-height:30px"):
        for i in range(len(cont)):
            if "--adapter" in cont[i]:
                with tag("pre", id="adapter_content"):
                    for j in range(i, i + 4):
                        text(cont[j] + "\n")
            elif "Top 5" in cont[i]:
                text(cont[i][4:] + "\n")
                with tag("div", style="line-height:20px"):
                    with tag("table", klass="customize", width="100%"):
                        for j in range(i + 1, i + 6):
                            linex = cont[j].split(" ")
                            with tag("tr"):
                                if j % 2 == 0:
                                    doc.attr(klass="odd")
                                else:
                                    doc.attr(klass="even")
                                for k in range(len(linex)):
                                    if (
                                        linex[k] != ""
                                        and ":" not in linex[k]
                                        and "=" not in linex[k]
                                    ):
                                        with tag("td", align="left"):
                                            text(linex[k])

    fin.close

def write_fastqc_report(doc, tag, text, line, report_dir, outputdir, max_sample=3):
    text(
        "The following is quality control file generated by FastQC. For more detailed information, please click the hyperlinks below."
    )
    sample_num = 0
    for root, dirs, files in os.walk(report_dir.getOutput("outputdir")):
        for report in files:
            if "1_fastqc.html" in report:
                sample_num += 1
                if (
                    sample_num > max_sample
                ):  # ignore the rest to shorten the report length
                    break
                with tag("div", id="fastqc_report_sub", klass="section level3"):
                    with tag("h3"):
                        text(
                            "Sample: "
                            + report.split("/")[-1].replace("_1_fastqc.html", "")
                        )
                    write_fastqc_report_contents(
                        doc,
                        tag,
                        text,
                        line,
                        os.path.join(report_dir.getOutput("outputdir"), report),
                        outputdir,
                    )
                    text(" ,  ")
                    write_fastqc_report_contents(
                        doc,
                        tag,
                        text,
                        line,
                        os.path.join(report_dir.getOutput("outputdir"), report.replace("1_fastqc.html", "2_fastqc.html")),
                        outputdir,
                    )

def write_fastqc_report_contents(doc, tag, text, line, report, outputdir):
    dstdir = os.path.join(outputdir, "./Fastq_Quality_Control/")
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    report_dir, report_name = os.path.split(report)
    shutil.copyfile(report, os.path.join(dstdir, report_name))
    with tag("a", href="Fastq_Quality_Control/" + report_name):
        text(report_name)
        
def write_deduplicate_report(doc, tag, text, line, report_dir, max_sample=3):
    sample_num = 0
    for report in report_dir.getOutput("reportOutput"):
        sample_num += 1
        if (
            sample_num > max_sample
        ):  # ignore the rest to shorten the report length
            break
        with tag("div", id="deduplicate_report_sub", klass="section level3"):
            with tag("h3"):
                text(
                    "Sample: "
                    + report.split("/")[-1].split(".")[0]
                )
            write_deduplicate_report_contents(
                doc,
                tag,
                text,
                line,
                report,
            )

def write_deduplicate_report_contents(doc, tag, text, line, report):
    f = open(report, "r")
    rs = f.readlines()
    with tag("div", style="line-height:30px"):
        with tag("ul", id="deduplicate_report_content"):
            line("li", rs[2])
            line("li", rs[5])
    f.close()

def write_rmduplicate_report(doc, tag, text, line, report_dir, max_sample=3):
    sample_num = 0
    with tag("div", id="rmduplicate_report_sub", klass="section level3"):
        with tag("div", style="line-height:20px"):
            with tag("table", klass="customize", width="100%"):
                with tag("tr"):
                    doc.attr(klass="odd")

                    with tag("td", align="left"):
                        text("SAMPLE_NAME")

                    with tag("td", align="left"):
                        text("READ_PAIRS_EXAMINED")

                    with tag("td", align="left"):
                        text("READ_PAIR_DUPLICATES")

                    with tag("td", align="left"):
                        text("PERCENT_DUPLICATION")

                    with tag("td", align="left"):
                        text("ESTIMATED_LIBRARY_SIZE")

                for report in report_dir.getOutput("metricsOutput"):
                    sample_num += 1
                    if (
                        sample_num > max_sample
                    ):  # ignore the rest to shorten the report length
                        break
                    with tag("tr"):
                        if sample_num % 2 == 0:
                            doc.attr(klass="odd")
                        else:
                            doc.attr(klass="even")
                        write_rmduplicate_report_contents(doc, tag, text, line, report)

def write_rmduplicate_report_contents(doc, tag, text, line, report):
    # read and locate rmdup.txt
    fin = open(report, "r", encoding="utf-8")
    conts = fin.readlines()
    cont = conts[7]
    data = cont.split("\t")

    # contents
    with tag("td", align="left"):
        text(report.split("/")[-1].replace("-rmdup.txt", ""))

    with tag("td", align="left"):
        text(data[2])

    with tag("td", align="left"):
        text(data[6])

    with tag("td", align="left"):
        text(data[8])

    with tag("td", align="left"):
        text(data[9])

    fin.close()

def write_fraglenplot_report(doc, tag, text, line, report_dir, outputdir, max_sample=3):
    sample_num = 0
    for report in report_dir.getOutput("plotOutput"):
        sample_num += 1
        if sample_num > max_sample:  # ignore the rest to shorten the report length
            break
        with tag("div", id="fraglenplot_report_sub", klass="section level3"):
            with tag("h3"):
                text("Sample: " + report.split("/")[-1].replace("_fraglen.png", ""))
            write_fraglenplot_report_contents(doc, tag, text, line, report, outputdir)

def write_fraglenplot_report_contents(doc, tag, text, line, report, outputdir):
    dstdir = outputdir + "/Fragment_Length/"
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    report_dir, report_name = os.path.split(report)
    dstfile = os.path.join(dstdir, report_name)
    shutil.copyfile(report, dstfile)
    doc.stag("img", src="Fragment_Length/" + report_name, alt=dstfile)

def write_OCF_report(doc, tag, text, line, report_dir, outputdir):
    report = report_dir.getOutput("plotOutput")
    with tag("div", id="OCF_report_sub", klass="section level3"):
        text("OCF value boxplot:")
    write_OCF_report_contents(doc, tag, text, line, report, outputdir)

def write_OCF_report_contents(doc, tag, text, line, report, outputdir):
    dstdir = outputdir + "/OCF/"
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    report_dir, report_name = os.path.split(report)
    dstfile = os.path.join(dstdir, report_name)
    shutil.copyfile(report, dstfile)
    doc.stag("img", src="OCF/" + report_name, alt=dstfile)
    
def write_CNV_report(doc, tag, text, line, report_dir, outputdir, label, max_sample=3):
    sample_num = 0
    casecorr = report_dir.getOutput("casereadplotOutput")
    ctrlcorr = report_dir.getOutput("ctrlreadplotOutput")
    with tag(
        "div",
        klass="corr",
    ):
        with tag("h2"):
            with tag("span", klass="section level2"):
                text("4.1 " + label[0] + " Read Count")
        for c in casecorr:
            sample_num += 1
            if sample_num > max_sample:  # ignore the rest to shorten the report length
                break
            with tag("div", id="case_corr_report_sub", klass="section level3"):
                with tag("h3"):
                    text("Sample: " + c.split("/")[-1].split(".")[0])
            write_CNV_report_contents(doc, tag, text, line, c, outputdir)
        sample_num = 0
        with tag("h2"):
            with tag("span", klass="section level2"):
                text("4.2 " + label[1] + " Read Count")
        for c in ctrlcorr:
            sample_num += 1
            if sample_num > max_sample:  # ignore the rest to shorten the report length
                break
            with tag("div", id="ctrl_corr_report_sub", klass="section level3"):
                with tag("h3"):
                    text("Sample: " + c.split("/")[-1].split(".")[0])
            write_CNV_report_contents(doc, tag, text, line, c, outputdir)
    with tag(
        "div",
        klass="z-core",
    ):
        with tag("h2"):
            with tag("span", klass="section level2"):
                text("4.3 Z-score plot")
        report = report_dir.getOutput("plotOutput")
        with tag("div", id="CNV_report_sub3", klass="section level3"):
            text("CNV z-score heatmap:")
        write_CNV_report_contents(doc, tag, text, line, report, outputdir)

def write_CNV_report_contents(doc, tag, text, line, report, outputdir):
    dstdir = outputdir + "/CNV/"
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    report_dir, report_name = os.path.split(report)
    dstfile = os.path.join(dstdir, report_name)
    shutil.copyfile(report, dstfile)
    doc.stag("img", src="CNV/" + report_name, alt=dstfile)