---
# Page settings
layout: default
keywords:
comments: true
image: paihia.jpg
# Hero section
title:  Chapter II
subtitle: >  
    A Predictive Theory of Allosteric Induction

# Author box
author:
    title: Summary
    title_url: ''
    external_url: false
    description: >
       TBD.
# Page navigation
page_nav:
    prev:
        content: Chapter 1
        url: chapter_01
    next:
        content: Chapter 3
        url: chapter_03
prefix: chapter_02
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_model
    - section_05_results
    - section_06_discussion
    - section_07_methods
---

**Published as ...**
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
