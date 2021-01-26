---
# Page settings
layout: default
keywords:
comments: true
image: alaska.jpg
# Hero section
title:  Chapter III
subtitle: >  
    First-principles prediction of the information processing capacity of a simple genetic circuit

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
        content: Chapter 2
        url: chapter_02
    next:
        content: Chapter 4
        url: chapter_04
prefix: chapter_03
contents:
    - section_01_header
    - section_02_abstract
    - section_03_introduction
    - section_04_minimal_model
    - section_05_distribution_moments
    - section_06_cell_cycle
    - section_07_maxent
    - section_08_channel_capacity
    - section_09_discussion
    - section_10_methods
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
