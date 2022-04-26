# Plotly

GlobalChem uses a lot of different graphing libraries and it can be cumbersome to write graphing templates over again. So to alleviate that problem, we provide a class option that will take in a figure and "beautify it" before you render to the output.&#x20;

The code is pretty simple but just pass your figure through the graphing template function and the template should be applied.&#x20;

**Apply the Template**

{% tabs %}
{% tab title="Code" %}
```
gce.apply_plotly_template(
            figure,
            x_title='X-Axis',
            y_title = 'Y-Axis',
            height = 500,
            width = 1000
)
```
{% endtab %}

{% tab title="Output" %}
![](<../.gitbook/assets/newplot (15).png>)
{% endtab %}
{% endtabs %}
