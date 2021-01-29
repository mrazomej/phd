---
# Page settings
layout: default
keywords:
comments: true
image: alaska.jpg
# Hero section
title:  Chapter IV
subtitle: >  
    Supporting Information for Tuning Transcriptional Regulation through Signaling: A Predictive Theory of Allosteric Induction


# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
        
# Page navigation
page_nav:
    prev:
        content: Chapter 3
        url: chapter_03
    next:
        content: Chapter 5
        url: chapter_05
prefix: chapter_04
contents:
    - section_01_header 
    - section_02_abstract
    - section_03_epsilonAI 
    - section_04_fugacity 
    - section_05_flow_gating 
    - section_06_microscopy 
    - section_07_sensitivity 
    - section_08_hill_fits 
    - section_09_global_fit 
    - section_10_oid 
    - section_11_comparison 
    - section_12_properties 
    - section_13_applications 
    - section_14_strainlist 
    - section_15_Nns_def 
    - section_16_steadystate
---

<hr/>
{% if page.contents %}
{% for val in page.contents %}
{% if jekyll.environment == production %}
{% include_relative {{site.doks.baseurl}}src/{{page.prefix}}/{{val}}.md %}
{% else %}
{% include_relative src/{{page.prefix}}/{{val}}.md %}
{% endif %}
{% endfor %}
{% endif %}

## References
