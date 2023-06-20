const Mathml2latex = require('mathml-to-latex');

const mathml = `
<math><mrow is="true"><mfrac is="true"><mrow is="true"><mi is="true">∂</mi><mi mathvariant="bold" is="true">c</mi></mrow><mrow is="true"><mi is="true">∂</mi><msup is="true"><mi is="true">t</mi><mo is="true">*</mo></msup></mrow></mfrac><mo linebreak="goodbreak" is="true">+</mo><mrow is="true"><mo stretchy="true" is="true">[</mo><mrow is="true"><mrow is="true"><mo stretchy="true" is="true">(</mo><mrow is="true"><mi mathvariant="bold" is="true">u</mi><mo linebreak="badbreak" is="true">−</mo><msub is="true"><mi mathvariant="bold" is="true">u</mi><mi is="true">c</mi></msub></mrow><mo stretchy="true" is="true">)</mo></mrow><mo is="true">·</mo><mi is="true">∇</mi></mrow><mo stretchy="true" is="true">]</mo></mrow><mi mathvariant="bold" is="true">c</mi><mo linebreak="goodbreak" is="true">=</mo><mrow is="true"><mo stretchy="true" is="true">(</mo><mrow is="true"><mi is="true">∇</mi><mi mathvariant="bold" is="true">u</mi></mrow><mo stretchy="true" is="true">)</mo></mrow><mo is="true">·</mo><mi mathvariant="bold" is="true">c</mi><mo linebreak="goodbreak" is="true">+</mo><mi mathvariant="bold" is="true">c</mi><msup is="true"><mrow is="true"><mo stretchy="true" is="true">(</mo><mrow is="true"><mi is="true">∇</mi><mi mathvariant="bold" is="true">u</mi></mrow><mo stretchy="true" is="true">)</mo></mrow><mi is="true">T</mi></msup><mo linebreak="goodbreak" is="true">−</mo><mfrac is="true"><msup is="true"><mrow is="true"><mrow is="true"><mi mathvariant="bold-italic" is="true">τ</mi></mrow></mrow><mi is="true">p</mi></msup><mrow is="true"><mi mathvariant="normal" is="true">W</mi><mi mathvariant="normal" is="true">e</mi></mrow></mfrac><mo linebreak="goodbreak" is="true">+</mo><mi is="true">κ</mi><mstyle mathvariant="normal" is="true"><mi is="true">Δ</mi></mstyle><mi mathvariant="bold" is="true">c</mi></mrow></math>
      `;


console.log(Mathml2latex.convert(mathml))