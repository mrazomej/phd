---
# Page settings
layout: default
keywords:
comments: true
image: paihia.jpg
# Hero section
title:  Chapter I
subtitle: >  
    Introduction

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
        TBD
# Page navigation
page_nav:
    prev:
        content: Abstract
        url: abstract
    next:
        content: Chapter 2
        url: chapter_02
prefix: chapter_01
contents:
    - section_00_title
    - section_01_introduction
    - section_02_gene_reg
    - section_03_information
---

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