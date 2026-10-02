import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}
const navTemp = {
  nav: [
    { text: 'Home', link: '/' },
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Background',
          items: [
            { text: 'Introduction', link: '/manual/introduction' },
            { text: 'Environments and the ecosystem', link: '/manual/ecosystem' },
            { text: 'Solving a heat balance', link: '/manual/heat_balance' },
            { text: 'Gradients, resistances and flows', link: '/manual/gradients' },
            { text: 'Units, dimensions and functional traits', link: '/manual/units_traits' },
          ]
        },
        { text: 'Solving',
          items: [
            { text: 'Temperature or metabolic rate', link: '/manual/solvers' },
            { text: 'Parameters', link: '/manual/parameters' },
          ]
        },
        { text: 'Processes',
          items: [
            { text: 'Radiation', link: '/manual/radiation' },
            { text: 'Convection and conduction', link: '/manual/convection_conduction' },
            { text: 'Evaporation and respiration', link: '/manual/evaporation_respiration' },
            { text: 'Metabolism', link: '/manual/metabolism' },
            { text: 'Insulation', link: '/manual/insulation' },
          ]
        },
        { text: 'Structure',
          items: [
            { text: 'Layers as a radial graph', link: '/manual/radial_layers' },
            { text: 'Bodies of many parts', link: '/manual/multipart' },
            { text: 'Differentiability and the NLP interface', link: '/manual/autodiff' },
          ]
        },
        { text: 'Reference',
          items: [
            { text: 'For NicheMapR users', link: '/manual/nichemapr' },
            { text: 'References', link: '/manual/references' },
          ]
        },
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'An ectotherm: body temperature', link: '/tutorials/ectotherm' },
        { text: 'An endotherm: metabolic rate', link: '/tutorials/endotherm' },
        { text: 'The endotherm, piece by piece', link: '/tutorials/components' },
        { text: 'A leaf', link: '/tutorials/leaf' },
        { text: 'The ellipsoid model', link: '/tutorials/ellipsoid' },
        { text: 'Back and belly: two halves', link: '/tutorials/two_parts' },
        { text: 'A human of many parts', link: '/tutorials/human' },
      ]
    },
    { text: 'API', link: '/api' }
  ],
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "Heat and water exchange between organisms and their environments, for biophysical ecology",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: false,
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Background', collapsed: false,
          items: [
            { text: 'Introduction', link: '/manual/introduction' },
            { text: 'Environments and the ecosystem', link: '/manual/ecosystem' },
            { text: 'Solving a heat balance', link: '/manual/heat_balance' },
            { text: 'Gradients, resistances and flows', link: '/manual/gradients' },
            { text: 'Units, dimensions and functional traits', link: '/manual/units_traits' },
          ]
        },
        { text: 'Solving', collapsed: false,
          items: [
            { text: 'Temperature or metabolic rate', link: '/manual/solvers' },
            { text: 'Parameters', link: '/manual/parameters' },
          ]
        },
        { text: 'Processes', collapsed: false,
          items: [
            { text: 'Radiation', link: '/manual/radiation' },
            { text: 'Convection and conduction', link: '/manual/convection_conduction' },
            { text: 'Evaporation and respiration', link: '/manual/evaporation_respiration' },
            { text: 'Metabolism', link: '/manual/metabolism' },
            { text: 'Insulation', link: '/manual/insulation' },
          ]
        },
        { text: 'Structure', collapsed: false,
          items: [
            { text: 'Layers as a radial graph', link: '/manual/radial_layers' },
            { text: 'Bodies of many parts', link: '/manual/multipart' },
            { text: 'Differentiability and the NLP interface', link: '/manual/autodiff' },
          ]
        },
        { text: 'Reference', collapsed: false,
          items: [
            { text: 'For NicheMapR users', link: '/manual/nichemapr' },
            { text: 'References', link: '/manual/references' },
          ]
        },
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'An ectotherm: body temperature', link: '/tutorials/ectotherm' },
        { text: 'An endotherm: metabolic rate', link: '/tutorials/endotherm' },
        { text: 'The endotherm, piece by piece', link: '/tutorials/components' },
        { text: 'A leaf', link: '/tutorials/leaf' },
        { text: 'The ellipsoid model', link: '/tutorials/ellipsoid' },
        { text: 'Back and belly: two halves', link: '/tutorials/two_parts' },
        { text: 'A human of many parts', link: '/tutorials/human' },
      ]
    },
    { text: 'API', link: '/api' }
    ],
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
